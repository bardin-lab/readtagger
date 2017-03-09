import argparse

import pysam
import six
import warnings

from cached_property import cached_property
from contextlib2 import ExitStack

from .cigar import alternative_alignment_cigar_is_better
from .bam_io import (
    BamAlignmentWriter as Writer,
    BamAlignmentReader as Reader,
)
from .tags import (
    BaseTag,
)

__VERSION__ = '0.3.6'


class SamTagProcessor(object):
    """Process SAM/AM file for tags of interest and keep a dict of readname, mate identity and tag in self.result."""

    def __init__(self, source_path, tag_prefix_self, tag_prefix_mate, tag_mate=True):
        """
        Process SAM/BAM at source path.

        :param source_path: Path of filesystem to SAM/BAM file
        :type source_path: basestring
        :param tag_prefix_self: Tag to use as indicating that the current read is aligned
        :type tag_prefix_self: basestring
        :param tag_prefix_mate: Tag to use as indicating that the current mate is aligned
        :type tag_prefix_mate: basestring
        :param tag_mate: Tag mate ?
        :type tag_mate: bool
        """
        self.tag_mate = tag_mate
        self.tag_prefix_self = tag_prefix_self
        self.tag_prefix_mate = tag_prefix_mate
        self.detail_tag_self = tag_prefix_self + 'D'
        self.detail_tag_mate = tag_prefix_mate + 'D'
        self.reference_tag_self = tag_prefix_self + 'R'
        self.reference_tag_mate = tag_prefix_mate + 'R'
        self.source_path = source_path
        self.tag_template = BaseTag(tid_to_reference_name=self.tid_to_reference_name)
        self.result = self.process_source()
        if self.tag_mate:
            self.add_mate()

    @cached_property
    def tid_to_reference_name(self):
        """Build dictionary for mapping tid to reference name."""
        with Reader(self.source_path, external_bin=None) as source:
            return {source.get_tid(rname): rname for rname in (d['SN'] for d in source.header['SQ'])}

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
        :return tags: namedtuple of information to add to tag
        :rtype namedtuple
        """
        return self.tag_template(tid=r.tid,
                                 reference_start=r.reference_start,
                                 cigar=r.cigar,
                                 is_reverse=r.is_reverse,
                                 mapq=r.mapping_quality,
                                 query_alignment_start=r.query_alignment_start,
                                 query_alignment_end=r.query_alignment_end)

    def is_taggable(self, r):
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
        ['readname'][Forward or Reverse][Self or Mate][Tag]
        """
        tag_d = {}
        with Reader(self.source_path) as source_alignment:
            for r in source_alignment:
                if self.is_taggable(r):
                    if r.query_name in tag_d:
                        tag_d[r.query_name][r.is_read1] = {'s': self.compute_tag(r)}
                    else:
                        tag_d[r.query_name] = {r.is_read1: {'s': self.compute_tag(r)}}
        return tag_d

    def add_mate(self):
        """Iterate over self.result and add mate information."""
        for tag_d in six.itervalues(self.result):
            for k, v in tag_d.items():
                if (not k) in tag_d:
                    v['m'] = tag_d[(not k)]['s']

    def format_tags(self, tags):
        """
        Format tags dictionary.

        >>> tags = {'m': (0, 362603, [(4, 2), (0, 37)], False, 4, 0, 37)}
        """
        formatted_tags = []
        for k, v in tags.items():
            if k == 's':
                formatted_tags.append((self.detail_tag_self, str(v)))
                formatted_tags.append((self.reference_tag_self, self.tid_to_reference_name[v[0]]))
            else:
                formatted_tags.append((self.detail_tag_mate, str(v)))
                formatted_tags.append(((self.reference_tag_mate, self.tid_to_reference_name[v.tid])))
        return formatted_tags

    def get_other_tag(self, other_r):
        """Take a read object `other_r` and fetch the annotation tag that has been processed in the instance."""
        other_r_reverse = other_r.is_read1
        tagged_mates = self.result.get(other_r.query_name)
        if tagged_mates:
            return tagged_mates.get(other_r_reverse, None)
        else:
            return None


class SamAnnotator(object):
    """Use list of SamTagProcessor instances to add tags to a BAM file."""

    def __init__(self,
                 annotate_file,
                 samtags,
                 output_path="test.bam",
                 allow_dovetailing=False,
                 discard_bad_alt_tag=True,
                 discard_if_proper_pair=False,
                 discarded_writer=None,
                 verified_writer=None):
        """
        Compare `samtags` with `annotate_file`.

        Produces a new alignment file at output_path.
        :param annotate_file: 'Path to bam/sam file'
        :type annotate_file: str
        :param samtags: list of SamTagProcessor instances
        :type samtags: List[SamTagProcessor]
        :param allow_dovetailing: Controls whether or not dovetailing should be allowed
        :type allow_dovetailing: bool
        """
        self.annotate_file = annotate_file
        self.output_path = output_path
        self.samtags = samtags
        if allow_dovetailing:
            self.max_proper_size = self.get_max_proper_pair_size(pysam.AlignmentFile(annotate_file))
            if not self.max_proper_size:
                allow_dovetailing = False
        self.process(allow_dovetailing, discard_bad_alt_tag, discard_if_proper_pair, discarded_writer, verified_writer)

    @cached_property
    def header(self):
        """Return SAM/BAM header."""
        with Reader(self.annotate_file, external_bin=None) as source:
            return source.header

    def process(self, allow_dovetailing=False, discard_bad_alt_tag=True, discard_if_proper_pair=False, discarded_writer=None, verified_writer=None):
        """Iterate over reads and fetch annotation from self.samtags."""
        kwds = {'reader': {Reader: {'path': self.annotate_file}},
                'main_writer': {Writer: {'path': self.output_path, 'header': self.header}},
                'discarded_writer': {Writer: {'path': discarded_writer, 'header': self.header}},
                'verified_writer': {Writer: {'path': verified_writer, 'header': self.header}}}
        args = {'allow_dovetailing': allow_dovetailing,
                'discard_bad_alt_tag': discard_bad_alt_tag,
                'discard_if_proper_pair': discard_if_proper_pair}
        with ExitStack() as stack:
            for arg, cls_arg_d in kwds.items():
                for cls, cls_args in cls_arg_d.items():
                    if cls_args['path']:
                        args[arg] = stack.enter_context(cls(**cls_args))
                    else:
                        args[arg] = None
            self._process(**args)

    def _process(self, reader, main_writer, discarded_writer, discard_if_proper_pair, verified_writer, allow_dovetailing, discard_bad_alt_tag):
        for read in reader:
            discarded_tags = []
            verified_tags = []
            verified_tag = None
            for samtag in self.samtags:
                alt_tag = samtag.get_other_tag(read)
                if alt_tag and discard_bad_alt_tag:
                    verified_tag = self.verify_alt_tag(read, alt_tag)
                    if discarded_writer and len(verified_tag) < len(alt_tag):
                        discarded_tag = samtag.format_tags({k: v for k, v in alt_tag.items() if k not in verified_tag})
                        discarded_tags.extend(discarded_tag)
                    alt_tag = verified_tag
                if alt_tag:
                    verified_tag = samtag.format_tags(alt_tag)
            if allow_dovetailing:
                read = self.allow_dovetailing(read, self.max_proper_size)
            if verified_tag and discard_if_proper_pair and read.is_proper_pair:
                discarded_tags.extend(verified_tag)
            elif verified_tag:
                verified_tags.extend(verified_tag)
            if discarded_tags:
                discarded_read = read.__copy__()
                discarded_read.tags += discarded_tags
                discarded_writer.write(discarded_read)
            if verified_tags:
                read.tags += verified_tags
                if verified_writer:
                    verified_writer.write(read)
            main_writer.write(read)

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
                if keep:
                    verified_alt_tag[s_or_m] = alt_tag
        return verified_alt_tag

    @classmethod
    def get_max_proper_pair_size(cls, alignment_file):
        """
        Iterate over the first 1000 properly paired records in alignment_file and get the maximum valid isize for a proper pair.

        :param alignment_file: pysam.AlignmentFile
        :type alignment_file: pysam.AlignmentFile
        :rtype int
        """
        isize = []
        for r in alignment_file:
            if r.is_proper_pair and not r.is_secondary and not r.is_supplementary:
                isize.append(abs(r.isize))
            if len(isize) == 1000:
                alignment_file.reset()
                return max(isize)
        alignment_file.reset()
        if isize:
            return max(isize)
        else:
            warnings.warn("Could not determine maximum allowed insert size for a proper pair. Are there any proper pairs in the input file?")
            return None

    @classmethod
    def allow_dovetailing(cls, read, max_proper_size=351, default_max_proper_size=351):
        """
        Manipulate is_proper_pair tag to allow dovetailing of reads.

        Precondition is read and mate have the same reference id, are within the maximum proper pair distance
        and are either in FR or RF orientation.
        :param read: aligned segment of pysam.AlignmentFile
        :type read: pysam.AlignedSegment
        :rtype pysam.AlignedSegment
        """
        if max_proper_size is None:
            warnings.warn("Using default maximum insert size of %d" % default_max_proper_size)
            max_proper_size = default_max_proper_size
        if not read.is_proper_pair and not read.is_reverse == read.mate_is_reverse and read.reference_id == read.mrnm and abs(read.isize) <= max_proper_size:
            read.is_proper_pair = True
        return read


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
    files_tags = zip(*parse_file_tags(args.annotate_with))
    samtags = [SamTagProcessor(filepath, tag_prefix_self=tag, tag_prefix_mate=tag_mate) for (filepath, tag, tag_mate) in files_tags]
    SamAnnotator(annotate_file=args.tag_file,
                 samtags=samtags,
                 output_path=args.output_file,
                 allow_dovetailing=args.allow_dovetailing,
                 discard_bad_alt_tag=not args.keep_suboptimal_alternate_tags,
                 discard_if_proper_pair=args.discard_if_proper_pair,
                 discarded_writer=args.write_discarded,
                 verified_writer=args.write_verified)


if __name__ == "__main__":
    main()
