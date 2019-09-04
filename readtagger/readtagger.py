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
    get_mean_read_length,
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
from tempfile import TemporaryDirectory

logger = logging.getLogger(__name__)
DEFAULT_CHUNK_SIZE = 10000
SELF = 's'
MATE = 'm'


class TagManager(object):
    """Orchestrates SamTagProcessor and SamAnnotator classes, holds reference to input and output files."""

    def __init__(self,
                 source_paths,
                 target_path,
                 output_path='test.bam',
                 cram=False,
                 discarded_path=None,
                 verified_path=None,
                 tag_mate=True,
                 allow_dovetailing=False,
                 max_proper_size=None,
                 discard_if_proper_pair=False,
                 discard_suboptimal_alternate_tags=True,
                 reference_fasta=None,
                 bwa_index=None,
                 tag_prefixes_self=('A',),
                 tag_prefixes_mate=('B',),
                 cores=1,
                 chunk_size='auto',
                 ):
        """Open input and output files and construct worker classes."""
        self.source_paths = source_paths
        self.source_paths_sorted = []
        self.annotate_path = target_path
        self.annotate_path_sorted = None
        self.output_path = output_path
        self.cram = cram
        self.discarded_path = discarded_path
        self.verified_path = verified_path
        self.tag_mate = tag_mate
        self.allow_dovetailing = allow_dovetailing
        self.max_proper_size = max_proper_size
        self.discard_if_proper_pair = discard_if_proper_pair
        self.discard_suboptimal_alternate_tags = discard_suboptimal_alternate_tags
        self.reference_fasta = reference_fasta
        self.bwa_index = bwa_index
        self.tag_prefixes_self = tag_prefixes_self
        self.tag_prefixes_mate = tag_prefixes_mate
        self.cores = cores
        self.chunk_size = chunk_size
        with TemporaryDirectory(prefix='TagManager_') as self.tempdir:
            self.setup_input_files()
            self.process()

    def setup_input_files(self):
        """Coordinate sort input files if necessary."""
        if is_file_coordinate_sorted(self.annotate_path):
            fd, path = tempfile.mkstemp()
            self.annotate_path_sorted = sort_bam(self.annotate_path, output=path, sort_order='queryname', threads=self.cores)
            os.close(fd)
        for source_path in self.source_paths:
            if is_file_coordinate_sorted(source_path):
                fd, path = tempfile.mkstemp()
                self.source_paths_sorted.append(sort_bam(source_path, output=path, sort_order='queryname', threads=self.cores))
                os.close(fd)
        self.annotate_path_sorted = self.annotate_path_sorted or self.annotate_path
        self.source_paths_sorted = self.source_paths_sorted or self.source_paths

    def setup_kwargs(self):
        """Set up arguments to pass to worker."""
        kwds = {}
        kwds['source_paths'] = self.source_paths_sorted
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
        kwds['tag_prefix_self'] = self.tag_prefixes_self
        kwds['tag_prefix_mate'] = self.tag_prefixes_mate
        kwds['source_headers'] = [pysam.AlignmentFile(p).header.to_dict() for p in self.source_paths_sorted]
        return kwds

    def setup_chunk_size(self):
        """Determine the chunk size."""
        if self.chunk_size == 'auto':
            # Adjust the chunk size based on read-length. We use the DEFAULT_CHUNK_SIZE for 200 nt reads
            # or less if the reads are longer (with a minimum of 10)
            mean_read_length = get_mean_read_length(self.annotate_path)
            chunk_normalization_factor = int(mean_read_length / 200)
            if chunk_normalization_factor == 0:
                chunk_normalization_factor = 1
            chunk_size = DEFAULT_CHUNK_SIZE / chunk_normalization_factor
            self.chunk_size = chunk_size if not chunk_size < 20 else 20
            logger.info("Chunk size is '%s', read size is '%s'", chunk_size, mean_read_length)

    def setup_source_file_splitting(self):
        """Find positions to which to jump in source path files."""
        kwds = self.setup_kwargs()
        mp_args = []
        pos_qname = get_queryname_positions(self.annotate_path_sorted, chunk_size=self.chunk_size, threads=self.cores)
        logger.info("Split annotate file %s into %d chunks" % (self.annotate_path_sorted, len(pos_qname)))
        last_qnames = [t[1] for t in pos_qname]
        starts_source = start_positions_for_last_qnames(self.source_paths_sorted[0], last_qnames=last_qnames, threads=self.cores)
        logger.info("Split source file %s into %d chunks" % (self.source_paths_sorted[0], len(pos_qname)))
        additional_source_starts = []
        if len(self.source_paths_sorted) > 1:
            for source_path in self.source_paths_sorted[1:]:
                additional_source_starts.append(start_positions_for_last_qnames(source_path, last_qnames=last_qnames, threads=self.cores))
        for i, ((start_annotate, qname), start_source) in enumerate(zip(pos_qname, starts_source)):
            args = kwds.copy()
            args['start_annotate'] = start_annotate
            start_source = [start_source]
            for sp in additional_source_starts:
                start_source.append(sp[i])
            args['start_source'] = start_source
            args['qname'] = qname
            args['chunk'] = i
            mp_args.append(args)

        return mp_args

    def process(self):
        """Create worker objects and stream pairs to SamAnnotator process method."""
        if self.allow_dovetailing and not self.max_proper_size:
            self.max_proper_size = get_max_proper_pair_size(self.annotate_path)

        if not self.bwa_index and self.reference_fasta:
            self.bwa_index, _ = make_bwa_index(reference_fasta=self.reference_fasta, dir=self.tempdir)

        self.setup_chunk_size()
        logger.info("Finding position at which to split input files")
        mp_args = self.setup_source_file_splitting()

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
    start_annotate = kwds['start_annotate']
    qname = kwds['qname']
    source_headers = kwds['source_headers']
    tempdir = kwds['tempdir']
    chunk = kwds['chunk']
    logger.info("Getting source reads for chunk %i", chunk)
    source_reads = [get_reads(p, start=s, last_qname=qname) for p, s in zip(kwds['source_paths'], kwds['start_source'])]
    logger.info("Getting target reads for chunk %i", chunk)
    annotate_reads = get_reads(kwds['annotate_path'], start=start_annotate, last_qname=qname)
    bwa_index = kwds.get('bwa_index')
    if bwa_index:
        TagSoftClip(source=annotate_reads, bwa_index=bwa_index, threads=2, min_clip_length=20)
    for reads in source_reads:
        AnnotateMateInformation(source=reads, target=annotate_reads)
    annotate_header = pysam.AlignmentFile(kwds['annotate_path']).header
    discarded_out = os.path.join(tempdir, "%s_discarded.bam" % chunk) if kwds['discarded_path'] else None
    discarded_writer = pysam.AlignmentFile(discarded_out, header=annotate_header, mode='wbu') if discarded_out else None
    verified_out = os.path.join(tempdir, "%s_verified.bam" % chunk) if kwds['verified_path'] else None
    verified_writer = pysam.AlignmentFile(verified_out, header=annotate_header, mode='wbu') if verified_out else None
    output_path = os.path.join(tempdir, "%s_output.bam" % chunk)
    output_writer = pysam.AlignmentFile(output_path, header=annotate_header, mode='wbu')
    samtag_p = [SamTagProcessor(source_bam=reads, header=headers, tag_mate=kwds['tag_mate']) for reads, headers in zip(source_reads, source_headers)]
    SamAnnotator(samtag_instances=samtag_p,
                 annotate_bam=annotate_reads,
                 output_writer=output_writer,
                 allow_dovetailing=kwds['allow_dovetailing'],
                 max_proper_size=kwds['max_proper_size'],
                 discard_suboptimal_alternate_tags=kwds['discard_suboptimal_alternate_tags'],
                 discard_if_proper_pair=kwds['discard_if_proper_pair'],
                 discarded_writer=discarded_writer,
                 verified_writer=verified_writer,
                 tag_prefixes_self=kwds['tag_prefix_self'],
                 tag_prefixes_mate=kwds['tag_prefix_mate'])
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
                    tag_d[r.query_name][r.is_read1] = {SELF: self.compute_tag(r)}
                else:
                    tag_d[r.query_name] = {r.is_read1: {SELF: self.compute_tag(r)}}
        return tag_d

    def add_mate(self):
        """Iterate over self.result and add mate information."""
        self._result = self.result
        self.result = {}
        for top_level_k, tag_d in self._result.items():
            new_tag_d = tag_d.copy()
            for k, v in tag_d.items():
                if (not k) in tag_d:
                    v[MATE] = tag_d[(not k)][SELF]
                    new_tag_d[k] = v
                else:
                    new_tag_d[(not k)] = {MATE: v[SELF]}
                self.result[top_level_k] = new_tag_d


class SamAnnotator(object):
    """Use list of SamTagProcessor instances to add tags to a BAM file."""

    def __init__(self,
                 samtag_instances,
                 annotate_bam,
                 output_writer=None,
                 allow_dovetailing=False,
                 max_proper_size=500,
                 discard_suboptimal_alternate_tags=True,
                 discard_if_proper_pair=False,
                 discarded_writer=None,
                 verified_writer=None,
                 tag_prefixes_self=('A',),
                 tag_prefixes_mate=('B',), ):
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
        self.samtag_instances = samtag_instances
        self.annotate_bam = annotate_bam
        self.output_writer = output_writer
        self.discarded_writer = discarded_writer
        self.verified_writer = verified_writer
        self.discard_suboptimal_alternate_tags = discard_suboptimal_alternate_tags
        self.discard_if_proper_pair = discard_if_proper_pair
        self.max_proper_size = max_proper_size
        self.allow_dovetailing = allow_dovetailing
        self.tag_prefixes_self = tag_prefixes_self
        self.tag_prefixes_mate = tag_prefixes_mate
        self.detail_tag_self = ["%sD" % t for t in self.tag_prefixes_self]
        self.detail_tag_mate = ["%sD" % t for t in self.tag_prefixes_mate]
        self.reference_tag_self = ["%sR" % t for t in self.tag_prefixes_self]
        self.reference_tag_mate = ["%sR" % t for t in self.tag_prefixes_mate]
        self.process()

    def process(self):
        """Process all reads in self.annotate_bam and self.samtag_instance."""
        for read in self.annotate_bam:
            if self.allow_dovetailing:
                read = allow_dovetailing(read, self.max_proper_size)
            discarded_tags = []
            verified_tags = []
            verified_tag = None
            for samtag_instance, detail_tag_self, details_tag_mate, reference_tag_self, reference_tag_mate in zip(self.samtag_instances,
                                                                                                                  self.detail_tag_self,
                                                                                                                  self.detail_tag_mate,
                                                                                                                  self.reference_tag_self,
                                                                                                                  self.reference_tag_mate):
                alt_tag = samtag_instance.result.get(read.query_name, {}).get(read.is_read1)
                if alt_tag and self.discard_suboptimal_alternate_tags:
                    # This is either not the correct read (unlikely because)
                    verified_tag = self.verify_alt_tag(read, alt_tag)
                    if self.discarded_writer and len(verified_tag) < len(alt_tag):
                        # we have more alt tags than verified tags,
                        # so we track the discarded tags
                        discarded_tag = self.format_tags({k: v for k, v in alt_tag.items() if k not in verified_tag},
                                                         detail_tag_self=detail_tag_self,
                                                         detail_tag_mate=details_tag_mate,
                                                         reference_tag_self=reference_tag_self,
                                                         reference_tag_mate=reference_tag_mate)
                        discarded_tags.extend(discarded_tag)
                    alt_tag = verified_tag
                if alt_tag:
                    verified_tag = self.format_tags(alt_tag,
                                                    detail_tag_self=detail_tag_self,
                                                    detail_tag_mate=details_tag_mate,
                                                    reference_tag_self=reference_tag_self,
                                                    reference_tag_mate=reference_tag_mate)
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

    def format_tags(self, tags, detail_tag_self, detail_tag_mate, reference_tag_self, reference_tag_mate):
        """
        Format tags dictionary.

        >>> tags = {'m': (0, 362603, [(4, 2), (0, 37)], False, 4, 0, 37)}
        """
        formatted_tags = []
        for k, v in tags.items():
            if k == 's':
                formatted_tags.append((detail_tag_self, str(v)))
                formatted_tags.append((reference_tag_self, v.reference_name()))
            else:
                formatted_tags.append((detail_tag_mate, str(v)))
                formatted_tags.append(((reference_tag_mate, v.reference_name())))
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
        for (alt_tag, self_or_mate, cigar) in self.get_tags_to_check(alt_tag=alt_tag, read=read):
            if cigar:
                # Can only check if read has cigar/alt_cigar
                same_orientation = alt_tag.is_reverse == (read.is_reverse if self_or_mate == SELF else read.mate_is_reverse)
                keep = alternative_alignment_cigar_is_better(current_cigar=cigar,
                                                             alternative_cigar=alt_tag.cigar,
                                                             same_orientation=same_orientation)
                if keep or (self_or_mate == MATE and not read.is_proper_pair):
                    # either the transposon tag is better, or we are on a different chromosome (perhaps a scaffold),
                    # in which case we still want to consider the alternative tag
                    verified_alt_tag[self_or_mate] = alt_tag
            elif self_or_mate == MATE:
                if read.mate_is_unmapped:
                    verified_alt_tag[self_or_mate] = alt_tag
        return verified_alt_tag

    def get_tags_to_check(self, alt_tag, read):
        """Determine the appropriate tag to check."""
        for self_or_mate, tag in alt_tag.items():
            if self_or_mate == SELF:
                yield (tag, SELF, read.cigar)
            else:
                if read.has_tag('MC'):
                    yield (tag, MATE, read.get_tag('MC'))
                else:
                    yield (tag, MATE, None)
