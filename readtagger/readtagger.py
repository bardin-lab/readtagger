import logging
import multiprocessing as mp
import os
import tempfile

import pysam

from .allow_dovetailing import (
    allow_dovetailing,
    get_max_proper_pair_size
)
from .bam_io import (
    BamAlignmentReader as Reader,
    get_queryname_positions,
    get_reads,
    is_file_coordinate_sorted,
    start_positions_for_last_qnames,
    merge_bam,
    sort_bam
)
from .bwa import make_bwa_index
from .cigar import alternative_alignment_cigar_is_better
from .mateoperations import AnnotateMateInformation
from .tags import (
    BaseTag,
    make_tag
)
from .tag_softclip import TagSoftClip

logger = logging.getLogger(__name__)
logging.basicConfig(format='%(asctime)s %(name)s %(levelname)s - %(message)s', level=logging.DEBUG)


class TagManager(object):
    """Orchestrates SamTagProcessor and SamAnnotator classes, holds reference to input and output files."""

    def __init__(self, source_path,
                 target_path,
                 output_path='test.bam',
                 discarded_path=None,
                 verified_path=None,
                 tag_mate=True,
                 allow_dovetailing=False,
                 max_proper_size=None,
                 discard_if_proper_pair=False,
                 discard_suboptimal_alternate_tags=True,
                 reference_fasta=None,
                 bwa_index=None,
                 tag_prefix_self='A',
                 tag_prefix_mate='B',
                 cores=1,
                 chunk_size=10000,
                 ):
        """Open input and output files and construct worker classes."""
        self.source_path = source_path
        self.source_path_sorted = None
        self.annotate_path = target_path
        self.annotate_path_sorted = None
        self.output_path = output_path
        self.discarded_path = discarded_path
        self.verified_path = verified_path
        self.tag_mate = tag_mate
        self.allow_dovetailing = allow_dovetailing
        self.max_proper_size = max_proper_size
        self.discard_if_proper_pair = discard_if_proper_pair
        self.discard_suboptimal_alternate_tags = discard_suboptimal_alternate_tags
        self.reference_fasta = reference_fasta
        self.bwa_index = bwa_index
        self.tag_prefix_self = tag_prefix_self
        self.tag_prefix_mate = tag_prefix_mate
        self.cores = cores
        self.chunk_size = chunk_size
        self.tempdir = tempfile.mkdtemp()
        self.setup_input_files()
        self.process()

    def setup_input_files(self):
        """Coordinate sort input files if necessary."""
        if is_file_coordinate_sorted(self.annotate_path):
            fd, path = tempfile.mkstemp()
            self.annotate_path_sorted = sort_bam(self.annotate_path, output=path, sort_order='queryname', threads=self.cores)
            os.close(fd)
        if is_file_coordinate_sorted(self.source_path):
            fd, path = tempfile.mkstemp()
            self.source_path_sorted = sort_bam(self.source_path, output=path, sort_order='queryname', threads=self.cores)
            os.close(fd)
        if not self.annotate_path_sorted:
            self.annotate_path_sorted = self.annotate_path
        if not self.source_path_sorted:
            self.source_path_sorted = self.source_path

    def process(self):
        """Create worker objects and stream pairs to SamAnnotator process method."""
        if self.allow_dovetailing and not self.max_proper_size:
            with Reader(self.annotate_path, external_bin=None) as source:
                self.max_proper_size = get_max_proper_pair_size(source)
        if not self.bwa_index and self.reference_fasta:
            tempdir = tempfile.mkdtemp()
            self.bwa_index, _ = make_bwa_index(reference_fasta=self.reference_fasta, dir=tempdir)
        kwds = {}
        kwds['source_path'] = self.source_path_sorted
        kwds['annotate_path'] = self.annotate_path_sorted
        kwds['discarded_path'] = self.discarded_path
        kwds['verified_path'] = self.verified_path
        kwds['bwa_index'] = self.bwa_index
        kwds['tempdir'] = self.tempdir
        kwds['tag_mate'] = self.tag_mate
        kwds['allow_dovetailing'] = self.allow_dovetailing
        kwds['max_proper_size'] = self.max_proper_size
        kwds['discard_suboptimal_alternate_tags'] = self.discard_suboptimal_alternate_tags
        kwds['discard_if_proper_pair'] = self.discard_if_proper_pair
        kwds['tag_prefix_self'] = self.tag_prefix_self
        kwds['tag_prefix_mate'] = self.tag_prefix_mate
        kwds['source_header'] = pysam.AlignmentFile(self.source_path_sorted).header
        logger.info("Finding position at which to split input files")
        pos_qname = get_queryname_positions(self.source_path_sorted, chunk_size=self.chunk_size)
        last_qnames = [t[1] for t in pos_qname]
        starts_annotate = start_positions_for_last_qnames(self.annotate_path_sorted, last_qnames=last_qnames)
        mp_args = []
        for (start_source, qname), start_annotate in zip(pos_qname, starts_annotate):
            args = kwds.copy()
            args['start_annotate'] = start_annotate
            args['start_source'] = start_source
            args['qname'] = qname
            mp_args.append(args)
        logger.info("Split processing into %d chunks", len(mp_args))
        if self.cores > 1:
            p = mp.Pool(self.cores)
            r = p.map_async(multiprocess_worker, mp_args)
            r = [o for o in r.get()]
        else:
            r = map(multiprocess_worker, mp_args)
        output_collection = []
        verified_collection = []
        discarded_collection = []
        logger.info("Writing output")
        for output, discarded, verified in r:
            output_collection.append(output)
            verified_collection.append(verified)
            discarded_collection.append(discarded)
        if any(discarded_collection):
            merge_bam(discarded_collection, template_bam=self.annotate_path, output_path=self.discarded_path)
            sort_bam(inpath=self.discarded_path, output=self.discarded_path, sort_order='coordinate', threads=self.cores)
        if any(verified_collection):
            merge_bam(verified_collection, template_bam=self.annotate_path, output_path=self.verified_path)
            sort_bam(inpath=self.verified_path, output=self.verified_path, sort_order='coordinate', threads=self.cores)
        merge_bam(output_collection, template_bam=self.annotate_path, output_path=self.output_path)
        sort_bam(self.output_path, output=self.output_path, sort_order='coordinate', threads=self.cores)
        logger.info("Processing finished")


def multiprocess_worker(kwds):
    """Process chunks of input bam files."""
    # source_bam can be a subset of annotate_bam
    start_source = kwds['start_source']
    start_annotate = kwds['start_annotate']
    qname = kwds['qname']
    source_header = kwds['source_header']
    tempdir = kwds['tempdir']
    source_reads = get_reads(kwds['source_path'], start=start_source, last_qname=qname)
    annotate_reads = get_reads(kwds['annotate_path'], start=start_annotate, last_qname=qname)
    bwa_index = kwds.get('bwa_index')
    if bwa_index:
        TagSoftClip(source=annotate_reads, bwa_index=bwa_index, threads=2, min_clip_length=20)
    AnnotateMateInformation(source=source_reads, target=annotate_reads)
    annotate_header = pysam.AlignmentFile(kwds['annotate_path']).header
    discarded_out = os.path.join(tempdir, "%s_discarded.bam" % qname) if kwds['discarded_path'] else None
    discarded_writer = pysam.AlignmentFile(discarded_out, header=annotate_header, mode='wbu') if discarded_out else None
    verified_out = os.path.join(tempdir, "%s_verified.bam" % qname) if kwds['verified_path'] else None
    verified_writer = pysam.AlignmentFile(verified_out, header=annotate_header, mode='wbu') if verified_out else None
    output_path = os.path.join(tempdir, "%s_output.bam" % qname)
    output_writer = pysam.AlignmentFile(output_path, header=annotate_header, mode='wbu')
    samtag_p = SamTagProcessor(source_bam=source_reads, header=source_header, tag_mate=kwds['tag_mate'])
    SamAnnotator(samtag_instance=samtag_p,
                 annotate_bam=annotate_reads,
                 output_writer=output_writer,
                 allow_dovetailing=kwds['allow_dovetailing'],
                 max_proper_size=kwds['max_proper_size'],
                 discard_suboptimal_alternate_tags=kwds['discard_suboptimal_alternate_tags'],
                 discard_if_proper_pair=kwds['discard_if_proper_pair'],
                 discarded_writer=discarded_writer,
                 verified_writer=verified_writer,
                 tag_prefix_self=kwds['tag_prefix_self'],
                 tag_prefix_mate=kwds['tag_prefix_mate'])
    if verified_writer:
        verified_writer.close()
    if discarded_writer:
        discarded_writer.close()
    output_writer.close()
    output_paths = [output_path, discarded_out, verified_out]
    return output_paths


class SamTagProcessor(object):
    """Process SAM/AM file for tags of interest and keep a dict of readname, mate identity and tag in self.result."""

    def __init__(self, source_bam, header, tag_mate=True):
        """
        Process SAM/BAM at source path.

        :param source_path: Path of filesystem to SAM/BAM file
        :type source_path: basestring

        :param tag_mate: Tag mate ?
        :type tag_mate: bool
        """
        self.tag_mate = tag_mate
        self.source_alignment = source_bam
        self.header = header
        self.template = BaseTag(header=self.header)
        self.result = self.process_source()
        if self.tag_mate:
            self.add_mate()

    def compute_tag(self, r):
        """
        Add tags for downstream processing.

        These are:
            - Reference ID
            - Reference start
            - CIGAR
            - Whether read us reverse
            - Mapping Quality
            - query_alignment_start
            - query_alignment_end

        :param r: AlignedSegment
        :rtype Tag
        """
        return make_tag(self.template, r=r)

    @staticmethod
    def is_taggable(r):
        """
        Decide if a read should be the source of a tag.

        :param r: AlignedSegment
        :type r: pysam.Alignedread
        """
        return not r.is_unmapped and not r.is_secondary and not r.is_supplementary and not r.is_qcfail

    def process_source(self):
        """
        Iterate over reads in alignment and construct dictionary.

        The dictionary struture is:
        ['readname'][Read1 or Read2][Self or Mate][Tag]
        """
        tag_d = {}
        for r in self.source_alignment:
            if self.is_taggable(r):
                if r.query_name in tag_d:
                    tag_d[r.query_name][r.is_read1] = {'s': self.compute_tag(r)}
                else:
                    tag_d[r.query_name] = {r.is_read1: {'s': self.compute_tag(r)}}
        return tag_d

    def add_mate(self):
        """Iterate over self.result and add mate information."""
        self._result = self.result
        self.result = {}
        for top_level_k, tag_d in self._result.items():
            new_tag_d = tag_d.copy()
            for k, v in tag_d.items():
                if (not k) in tag_d:
                    v['m'] = tag_d[(not k)]['s']
                    new_tag_d[k] = v
                else:
                    new_tag_d[(not k)] = {'m': v['s']}
                self.result[top_level_k] = new_tag_d


class SamAnnotator(object):
    """Use list of SamTagProcessor instances to add tags to a BAM file."""

    def __init__(self,
                 samtag_instance,
                 annotate_bam,
                 output_writer=None,
                 allow_dovetailing=False,
                 max_proper_size=500,
                 discard_suboptimal_alternate_tags=True,
                 discard_if_proper_pair=False,
                 discarded_writer=None,
                 verified_writer=None,
                 tag_prefix_self='A',
                 tag_prefix_mate='B', ):
        """
        Compare `samtags` with `annotate_file`.

        Produces a new alignment file at output_path.
        :param annotate_file: pysam.AlignmentFile for bam file to annotate
        :type annotate_file: pysam.AlignmentFile
        :param output_writer: pysam.AlignmentFile for output bam file
        :type output_writer: pysam.AlignmentFile
        :param discarded_writer: pysam.AlignmentFile for bam file containing reads with discarded tags
        :type discarded_writer: pysam.AlignmentFile
        :param verified_writer: pysam.AlignmentFile for output bam file containing only reads with verified tags
        :type verified_writer: pysam.AlignmentFile
        :param samtags: list of SamTagProcessor instances
        :type samtags: List[SamTagProcessor]
        :param allow_dovetailing: Controls whether or not dovetailing should be allowed
        :type allow_dovetailing: bool
        :param tag_prefix_self: Tag to use as indicating that the current read is aligned
        :type tag_prefix_self: basestring
        :param tag_prefix_mate: Tag to use as indicating that the current mate is aligned
        :type tag_prefix_mate: basestring
        """
        self.samtag_instance = samtag_instance
        self.annotate_bam = annotate_bam
        self.output_writer = output_writer
        self.discarded_writer = discarded_writer
        self.verified_writer = verified_writer
        self.discard_suboptimal_alternate_tags = discard_suboptimal_alternate_tags
        self.discard_if_proper_pair = discard_if_proper_pair
        self.max_proper_size = max_proper_size
        self.allow_dovetailing = allow_dovetailing
        self.tag_prefix_self = tag_prefix_self
        self.tag_prefix_mate = tag_prefix_mate
        self.detail_tag_self = tag_prefix_self + 'D'
        self.detail_tag_mate = tag_prefix_mate + 'D'
        self.reference_tag_self = tag_prefix_self + 'R'
        self.reference_tag_mate = tag_prefix_mate + 'R'
        self.process()

    def process(self):
        """Process all reads in self.annotate_bam and self.samtag_instance."""
        for read in self.annotate_bam:
            if self.allow_dovetailing:
                read = allow_dovetailing(read, self.max_proper_size)
            discarded_tags = []
            verified_tags = []
            verified_tag = None
            alt_tag = self.samtag_instance.result.get(read.query_name, {}).get(read.is_read1)
            if alt_tag and self.discard_suboptimal_alternate_tags:
                # This is either not the correct read (unlikely because)
                verified_tag = self.verify_alt_tag(read, alt_tag)
                if self.discarded_writer and len(verified_tag) < len(alt_tag):
                    # we have more alt tags than verified tags,
                    # so we track the discarded tags
                    discarded_tag = self.format_tags({k: v for k, v in alt_tag.items() if k not in verified_tag})
                    discarded_tags.extend(discarded_tag)
                alt_tag = verified_tag
            if alt_tag:
                verified_tag = self.format_tags(alt_tag)
            if verified_tag and self.discard_if_proper_pair and read.is_proper_pair:
                discarded_tags.extend(verified_tag)
            elif verified_tag:
                verified_tags.extend(verified_tag)
            if discarded_tags:
                discarded_read = read.__copy__()
                discarded_read.tags += discarded_tags
                if self.discarded_writer:
                    self.discarded_writer.write(discarded_read)
            if verified_tags:
                read.tags += verified_tags
                if self.verified_writer:
                    self.verified_writer.write(read)
            self.output_writer.write(read)

    def format_tags(self, tags):
        """
        Format tags dictionary.

        >>> tags = {'m': (0, 362603, [(4, 2), (0, 37)], False, 4, 0, 37)}
        """
        formatted_tags = []
        for k, v in tags.items():
            if k == 's':
                formatted_tags.append((self.detail_tag_self, str(v)))
                formatted_tags.append((self.reference_tag_self, v.reference_name()))
            else:
                formatted_tags.append((self.detail_tag_mate, str(v)))
                formatted_tags.append(((self.reference_tag_mate, v.reference_name())))
        return formatted_tags

    def verify_alt_tag(self, read, alt_tag):
        """
        Take a read and verify that the alternative tags are really better.

        :param read: A pysam read
        :type read: pysam.AlignedRead
        :return: read
        :type read: pysam.AlignedRead
        """
        # TODO: Make this available as a standalone function by explcitly passing in the tags to look at
        verified_alt_tag = {}
        if read.is_unmapped:
            return verified_alt_tag
        tags_to_check = ((v, 's', read.cigar) if k == 's' else (v, 'm', read.get_tag('MC')) if read.has_tag('MC') else (None, None, None) for k, v in
                         alt_tag.items())
        for (alt_tag, s_or_m, cigar) in tags_to_check:
            if cigar:
                # Can only check if read has cigar/alt_cigar
                same_orientation = alt_tag.is_reverse == (read.is_reverse if s_or_m == 's' else read.mate_is_reverse)
                keep = alternative_alignment_cigar_is_better(current_cigar=cigar,
                                                             alternative_cigar=alt_tag.cigar,
                                                             same_orientation=same_orientation)
                if keep or (s_or_m == 'm' and not read.is_proper_pair):
                    # either the transposon tag is better, or we are on a different chromosome (perhaps a scaffold),
                    # in which case we still want to consider the alternative tag
                    verified_alt_tag[s_or_m] = alt_tag
        return verified_alt_tag
