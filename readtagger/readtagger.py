import argparse

from contextlib2 import ExitStack

from .allow_dovetailing import (
    allow_dovetailing,
    get_max_proper_pair_size
)
from .cigar import alternative_alignment_cigar_is_better
from .bam_io import (
    BamAlignmentWriter as Writer,
    BamAlignmentReader as Reader,
)
from .tags import Tag

__VERSION__ = '0.3.14'


class TagManager(object):
    """Orchestrates SamTagProcessor and SamAnnotator classes, holds reference to input and output files."""

    def __init__(self, source_path,
                 annotate_path,
                 output_path='test.bam',
                 discarded_path=None,
                 verified_path=None,
                 tag_mate=True,
                 allow_dovetailing=False,
                 max_proper_size=None,
                 discard_if_proper_pair=False,
                 discard_bad_alt_tag=True,
                 tag_prefix_self='A',
                 tag_prefix_mate='B',
                 ):
        """Open input and output files and construct worker classes."""
        self.source_path = source_path
        self.annotate_path = annotate_path
        self.output_path = output_path
        self.discarded_path = discarded_path
        self.verified_path = verified_path
        self.tag_mate = tag_mate
        self.allow_dovetailing = allow_dovetailing
        self.max_proper_size = max_proper_size
        self.discard_if_proper_pair = discard_if_proper_pair
        self.discard_bad_alt_tag = discard_bad_alt_tag
        self.tag_prefix_self = tag_prefix_self
        self.tag_prefix_mate = tag_prefix_mate
        self.process()

    def process(self):
        """Create worker objects and stream pairs to SamAnnotator process method."""
        if self.allow_dovetailing and not self.max_proper_size:
            with Reader(self.source_path, external_bin=None) as source:
                self.max_proper_size = get_max_proper_pair_size(source)
        with Reader(self.annotate_path) as annotate:
            annotate_header = annotate.header
        kwds = {'source_bam': {Reader: {'path': self.source_path, 'sort_order': 'queryname'}},
                'annotate_bam': {Reader: {'path': self.annotate_path, 'sort_order': 'queryname'}},
                'output_writer': {Writer: {'path': self.output_path, 'header': annotate_header}},
                'discarded_writer': {Writer: {'path': self.discarded_path, 'header': annotate_header}},
                'verified_writer': {Writer: {'path': self.verified_path, 'header': annotate_header}}}
        args = {}
        with ExitStack() as stack:
            for arg, cls_arg_d in kwds.items():
                for cls, cls_args in cls_arg_d.items():
                    if cls_args['path']:
                        args[arg] = stack.enter_context(cls(**cls_args))
                    else:
                        args[arg] = None
            samtag_p = SamTagProcessor(source_bam=args['source_bam'], tag_mate=self.tag_mate)
            sam_annotator = SamAnnotator(annotate_bam=args['annotate_bam'],
                                         output_writer=args['output_writer'],
                                         allow_dovetailing=self.allow_dovetailing,
                                         max_proper_size=self.max_proper_size,
                                         discard_bad_alt_tag=self.discard_bad_alt_tag,
                                         discard_if_proper_pair=self.discard_if_proper_pair,
                                         discarded_writer=args['discarded_writer'],
                                         verified_writer=args['verified_writer'],
                                         tag_prefix_self=self.tag_prefix_self,
                                         tag_prefix_mate=self.tag_prefix_mate
                                         )
            for qname, tag_d in samtag_p.process():
                sam_annotator.process(qname=qname, tag_d=tag_d)


class SamTagProcessor(object):
    """Process SAM/AM file for tags of interest and keep a dict of readname, mate identity and tag in self.result."""

    def __init__(self, source_bam, tag_mate=True):
        """
        Process SAM/BAM at source path.

        :param source_path: Path of filesystem to SAM/BAM file
        :type source_path: basestring

        :param tag_mate: Tag mate ?
        :type tag_mate: bool
        """
        self.tag_mate = tag_mate
        self.source_alignment = source_bam
        self.header = source_bam.header

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
        return Tag.from_read(r, header=self.header)

    @staticmethod
    def is_taggable(r):
        """
        Decide if a read should be the source of a tag.

        :param r: AlignedSegment
        :type r: pysam.Alignedread
        """
        return not r.is_unmapped and not r.is_secondary and not r.is_supplementary and not r.is_qcfail

    def process(self):
        """Generator that will return qname and tag_d."""
        reads = []
        qname = ''
        for read in self.source_alignment:
            if self.is_taggable(read):
                if qname != read.query_name and reads and qname:
                    # We got a new read, so we process the revious reads
                    tag_d = self.process_current_reads(reads)
                    if self.tag_mate:
                        tag_d = self.add_mate(tag_d)
                    yield qname, tag_d
                    reads = []
                reads.append(read)
            qname = read.query_name
        if reads:
            tag_d = self.process_current_reads(reads)
            if self.tag_mate:
                tag_d = self.add_mate(tag_d)
            yield qname, tag_d

    def process_current_reads(self, current_reads):
        """
        Iterate over reads in alignment and construct dictionary.

        The dictionary struture is:
        ['readname'][Forward or Reverse][Self or Mate][Tag]
        """
        tag_d = {}
        for r in current_reads:
            if self.is_taggable(r):
                tag_d[r.is_read1] = {'s': self.compute_tag(r)}
        return tag_d

    def add_mate(self, tag_d):
        """Iterate over tag_d and add mate information."""
        new_tag_d = tag_d.copy()
        for k, v in tag_d.items():
            if (not k) in tag_d:
                v['m'] = tag_d[(not k)]['s']
                new_tag_d[k] = v
            else:
                new_tag_d[(not k)] = {'m': v['s']}
        return new_tag_d


class SamAnnotator(object):
    """Use list of SamTagProcessor instances to add tags to a BAM file."""

    def __init__(self,
                 annotate_bam,
                 output_writer=None,
                 allow_dovetailing=False,
                 max_proper_size=500,
                 discard_bad_alt_tag=True,
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
        self.annotate_bam = annotate_bam
        self.header = annotate_bam.header
        self.output_writer = output_writer
        self.discarded_writer = discarded_writer
        self.verified_writer = verified_writer
        self.discard_bad_alt_tag = discard_bad_alt_tag
        self.discard_if_proper_pair = discard_if_proper_pair
        self.max_proper_size = max_proper_size
        self.allow_dovetailing = allow_dovetailing
        self.tag_prefix_self = tag_prefix_self
        self.tag_prefix_mate = tag_prefix_mate
        self.detail_tag_self = tag_prefix_self + 'D'
        self.detail_tag_mate = tag_prefix_mate + 'D'
        self.reference_tag_self = tag_prefix_self + 'R'
        self.reference_tag_mate = tag_prefix_mate + 'R'
        self._read = None

    def process(self, qname, tag_d):
        """Iterate through reads in self.annotate_bam and process those reads that have the same name as qname."""
        reads = []
        if self._read:
            # We saved the last read from the previous process call
            if self._read.query_name == qname:
                reads.append(self._read)
        for read in self.annotate_bam:
            if read.query_name == qname:
                reads.append(read)
            else:
                if not reads:
                    # That should be the first call, when no reads have been found yet and
                    # we haven't reached the correct reads yet
                    continue
                else:
                    # We read past the last read with the same name
                    return self._process(reads, tag_d)
        # We should reach this only at the end of a BAM file in rare instances
        return self._process(reads, tag_d)

    def _process(self, reads, tag_d):
        for read in reads:
            if allow_dovetailing:
                read = allow_dovetailing(read, self.max_proper_size)
            discarded_tags = []
            verified_tags = []
            verified_tag = None
            alt_tag = tag_d.get(read.is_reverse)
            if alt_tag and self.discard_bad_alt_tag:
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
                formatted_tags.append((self.detail_tag_self, v.to_string()))
                formatted_tags.append((self.reference_tag_self, v.reference_name))
            else:
                formatted_tags.append((self.detail_tag_mate, v.to_string()))
                formatted_tags.append(((self.reference_tag_mate, v.reference_name)))
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


def parse_file_tags(filetags):
    """
    Parse list of filetags from commandline.

    :param filetags: list of strings with filepath.
                     optionally appended by the first letter that should be used for read and mate
    :return: annotate_with, tag_prefix, tag_prefix_mate

    >>> filetags = ['file_a:A:B', 'file_b:C:D', 'file_c']
    >>> annotate_with, tag_prefix, tag_prefix_mate = parse_file_tags(filetags)
    >>> annotate_with == ['file_a', 'file_b', 'file_c'] and tag_prefix == ['A', 'C', 'A'] and tag_prefix_mate == ['B', 'D', 'B']
    True
    >>>
    """
    annotate_with = []
    tag_prefix = []
    tag_prefix_mate = []
    if not isinstance(filetags, list):
        filetags = [filetags]
    for filetag in filetags:
        if ':' in filetag:
            filepath, tag, tag_mate = filetag.split(':')
            annotate_with.append(filepath)
            tag_prefix.append(tag.upper())
            tag_prefix_mate.append(tag_mate.upper())
        else:
            annotate_with.append(filetag)
            tag_prefix.append('A')  # Default is A for read, B for mate
            tag_prefix_mate.append('B')
    return annotate_with, tag_prefix, tag_prefix_mate


def parse_args():
    """Parse commandline arguments."""
    p = argparse.ArgumentParser(description="Tag reads in an alignment file based on other alignment files",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-t', '--tag_file', help="Tag reads in this file.", required=True)
    p.add_argument('-a', '--annotate_with',
                   help="Tag reads in readfile if reads are aligned in these files."
                        "Append `:A:B` to tag first letter of tag describing read as A, "
                        "and first letter of tag describing the mate as B",
                   nargs="+",
                   required=True)
    p.add_argument('-o', '--output_file', help="Write bam file to this path", required=True)
    p.add_argument('-d', '--allow_dovetailing',
                   action='store_true',
                   help="Sets the proper pair flag (0x0002) to true if reads dovetail [reads reach into or surpass the mate sequence].")
    p.add_argument('-dp', '--discard_if_proper_pair', action='store_true', help="Discard an alternative flag if the current read is in a proper pair.")
    p.add_argument('-k', '--keep_suboptimal_alternate_tags', action='store_true',
                   help="By default cigarstrings of the alternative tags are compared and alternates that are not explaining the current cigar "
                        "strings are discarded. Use this option to keep the alternative tags (effectively restoring the behaviour of readtagger < 0.1.4)")
    p.add_argument('-wd', '--write_discarded', default=False, required=False, help="Write discarded reads into separate file")
    p.add_argument('-wv', '--write_verified', default=False, required=False,
                   help="Write verified reads into separate file")
    p.add_argument('--version', action='version', version=__VERSION__)
    return p.parse_args()


def main(args=None):
    """Main entrypoint."""
    if not args:
        args = parse_args()
    filepath, tag, tag_mate = list(zip(*parse_file_tags(args.annotate_with)))[0]
    TagManager(
        source_path=filepath,
        annotate_path=args.tag_file,
        output_path=args.output_file,
        discarded_path=args.write_discarded,
        verified_path=args.write_verified,
        tag_mate=True,
        allow_dovetailing=args.allow_dovetailing,
        discard_if_proper_pair=args.discard_if_proper_pair,
        discard_bad_alt_tag=not args.keep_suboptimal_alternate_tags,
        tag_prefix_self=tag,
        tag_prefix_mate=tag_mate)


if __name__ == "__main__":
    main()
