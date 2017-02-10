import argparse
import parse
import pysam
import six

from .cigar import alternative_alignment_cigar_is_better


TAG_TEMPLATE = "R:{ref},POS:{pos:d},QSTART:{qstart:d},QEND:{qend:d},CIGAR:{cigar},S:{sense},MQ:{mq:d}"
TAG_TEMPLATE_COMPILED = parse.compile(TAG_TEMPLATE)


class SamTagProcessor(object):

    tag_template = TAG_TEMPLATE

    def __init__(self, source_path, tag_prefix_self, tag_prefix_mate, tag_mate=True):
        self.tag_mate = tag_mate
        self.tag_prefix_self = tag_prefix_self
        self.tag_prefix_mate = tag_prefix_mate
        self.detail_tag_self = tag_prefix_self + 'D'
        self.detail_tag_mate = tag_prefix_mate + 'D'
        self.reference_tag_self = tag_prefix_self + 'R'
        self.reference_tag_mate = tag_prefix_mate + 'R'
        self.source_alignment = pysam.AlignmentFile(source_path)
        self.result = self.process_source()
        if self.tag_mate:
            self.add_mate()

    def compute_tag(self, r):
        """
        Adds tags for downstream processing:
            - Reference
            - Pos
            - Sense
            - Aligned position
            - MD
        Returns a tag of the form:
        [('AA'), ('R:gypsy,POS:1,MD:107S35M,S:AS')]
        'AA' stands for alternative alignment,
        'MA' stands for mate alignment.
        :param r: AlignedSegment
        """
        tags = (self.source_alignment.get_reference_name(r.tid),
                r.reference_start,
                r.cigarstring,
                'AS' if r.is_reverse else 'S',
                r.mapping_quality,
                r.qstart,
                r.qend)
        return tags

    def is_taggable(self, r):
        """
        Decide if a read should be tagged
        :param r: AlignedSegment
        :type r: pysam.Alignedread
        """
        return not r.is_qcfail and not r.is_secondary and not r.is_supplementary and not r.is_unmapped

    def process_source(self):
        tag_d = {}
        for r in self.source_alignment:
            if self.is_taggable(r):
                fw_rev = 'r1' if r.is_read1 else 'r2'
                if r.query_name in tag_d:
                    tag_d[r.query_name][fw_rev] = {'s': self.compute_tag(r)}
                else:
                    tag_d[r.query_name] = {fw_rev: {'s': self.compute_tag(r)}}
        return tag_d

    def add_mate(self):
        for qname, tag_d in six.iteritems(self.result):
            if len(tag_d) == 2:
                # Both mates aligned, add mate tag
                tag_d['r1']['m'] = tag_d['r2']['s']
                tag_d['r2']['m'] = tag_d['r1']['s']
            elif len(tag_d) == 1:
                # Only one of the mates mapped, so we fill the mate with its mate tag
                if 'r1' in tag_d:
                    tag_d['r2'] = {'m': tag_d['r1']['s']}
                else:
                    tag_d['r1'] = {'m': tag_d['r2']['s']}
            else:
                continue
            self.result[qname] = tag_d

    def format_tag_value(self, t):
        return self.tag_template.format(ref=t[0], pos=t[1], cigar=t[2], sense=t[3], mq=t[4], qstart=t[5], qend=t[6])

    def format_tags(self, tags):
        formatted_tags = []
        if 'm' in tags:
            formatted_tags.append((self.detail_tag_mate, self.format_tag_value(tags['m'])))
            formatted_tags.append((self.reference_tag_mate, tags['m'][0]))
        if 's' in tags:
            formatted_tags.append((self.detail_tag_self, self.format_tag_value(tags['s'])))
            formatted_tags.append((self.reference_tag_self, tags['s'][0]))
        return formatted_tags

    def get_other_tag(self, other_r):
        """
        convenience method that takes a read object `other_r` and fetches the
        annotation tag that has been processed in the instance
        """
        other_r_reverse = 'r1' if other_r.is_read1 else 'r2'
        tagged_mates = self.result.get(other_r.query_name)
        if tagged_mates:
            return self.format_tags(tagged_mates[other_r_reverse])
        else:
            return None

class SamAnnotator(object):

    tag_template_compiled = TAG_TEMPLATE_COMPILED

    def __init__(self, annotate_file, samtags, output_path="test.bam", allow_dovetailing=False, discard_bad_alt_tag=True):
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
        self.annotate_file = pysam.AlignmentFile(annotate_file)
        self.output = pysam.AlignmentFile(output_path, 'wb', template=self.annotate_file)
        self.samtags = samtags
        self.detail_tag_self = self.samtags[0].detail_tag_self  # urgs ... probably bad design here :(.
        self.detail_tag_mate = self.samtags[0].detail_tag_mate
        self.reference_tag_self = self.samtags[0].reference_tag_self
        self.reference_tag_mate = self.samtags[0].reference_tag_mate
        self.process(allow_dovetailing, discard_bad_alt_tag)

    def process(self, allow_dovetailing=False, discard_bad_alt_tag=True):
        if allow_dovetailing:
            max_proper_size = self.get_max_proper_pair_size(self.annotate_file)
            if not max_proper_size:
                allow_dovetailing = False
        for read in self.annotate_file:
            for samtag in self.samtags:
                annotated_tag = samtag.get_other_tag(read)
                if annotated_tag:
                    read.tags = annotated_tag + read.tags
                    if discard_bad_alt_tag:
                        read = self.verify_alt_tag(read)
            if allow_dovetailing:
                read = self.allow_dovetailing(read, max_proper_size)
            self.output.write(read)
        self.output.close()

    def verify_alt_tag(self, read):
        """
        Take a read and verify that the alternative tags are really better
        :param read: A pysam read
        :type read: pysam.AlignedRead
        :return: read
        :type read: pysam.AlignedRead
        """
        # TODO: Make this available as a standalone function by explcitly passing in the tags to look at
        tags_to_check = []
        if read.has_tag(self.detail_tag_self):
            tags_to_check.append((self.detail_tag_self, self.reference_tag_self, read.cigarstring))
        if read.has_tag(self.detail_tag_mate) and read.has_tag('MC'):
            mc_tag = read.get_tag('MC')
            tags_to_check.append((self.detail_tag_mate, self.reference_tag_mate, mc_tag))
        for (alt_tag, alt_r_tag, cigar) in tags_to_check:
            tag = self.tag_template_compiled.parse(read.get_tag(alt_tag))
            alt_is_reverse = True if tag['sense'] == 'AS' else False
            keep = alternative_alignment_cigar_is_better(current_cigar=cigar,
                                                         alternative_cigar=tag['cigar'],
                                                         same_orientation=alt_is_reverse == read.is_reverse)
            if not keep:
                read.set_tag(tag=alt_tag, value=None)
                read.set_tag(tag=alt_r_tag, value=None)
        return read


    @classmethod
    def get_max_proper_pair_size(cls, alignment_file):
        """
        iterate over the first 1000 properly paired records in alignment_file
        and get the maximum valid isize for a proper pair.
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
            return None

    @classmethod
    def allow_dovetailing(cls, read, max_proper_size=351):
        """
        Manipulates is_proper_pair tag to allow dovetailing of reads.
        Precondition is read and mate have the same reference id, are within the maximum proper pair distance
        and are either in FR or RF orientation.
        :param read: aligned segment of pysam.AlignmentFile
        :type read: pysam.AlignedSegment
        :rtype pysam.AlignedSegment
        """
        if not read.is_proper_pair and not read.is_reverse == read.mate_is_reverse and read.reference_id == read.mrnm and abs(read.isize) <= max_proper_size:
            read.is_proper_pair = True
        return read


def parse_file_tags(filetags):
    """
    :param filetags: string with filepath.
                     optionally appended by the first letter that should be used for read and mate
    :return: annotate_with, tag_prefix, tag_prefix_mate

    >>> filetags = ['file_a:A:B', 'file_b:C:D', 'file_c']
    >>> annotate_with, tag_prefix, tag_prefix_mate = parse_file_tags(filetags)
    >>> annotate_with == ['file_a', 'file_b', 'file_c'] and tag_prefix == ['A', 'C', 'R'] and tag_prefix_mate == ['B', 'D', 'M']
    True
    >>>
    """
    annotate_with = []
    tag_prefix = []
    tag_prefix_mate = []
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
    p = argparse.ArgumentParser(description="Tag reads in an alignment file based on other alignment files",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-t', '--tag_file', help="Tag reads in this file.", required=True)
    p.add_argument('-a', '--annotate_with',
                   help="Tag reads in readfile if reads are aligned in these files."
                        "Append `:A:B` to tag first letter of tag describing read as A, "
                        "and first letter of tag describing the mate as B",
                   nargs = "+",
                   required=True)
    p.add_argument('-o', '--output_file', help="Write bam file to this path", required=True)
    p.add_argument('-d', '--allow_dovetailing',
                   action='store_true',
                   help="Sets the proper pair flag (0x0002) to true if reads dovetail [reads reach into or surpass the mate sequence].")
    p.add_argument('-k', '--keep_suboptimal_alternate_tags', action='store_true',
                   help="By default cigar strings of the alternative tags are compared and alternates that are not explaining the current cigar strings are discarded."
                        "Use this option to keep the alternative tags (effectively restoring the behaviour of tag_reads < 0.1.4)")
    return p.parse_args()

def main():
    args = parse_args()
    files_tags = zip(*parse_file_tags(args.annotate_with))
    samtags = [SamTagProcessor(filepath, tag_prefix_self=tag, tag_prefix_mate=tag_mate) for (filepath, tag, tag_mate) in files_tags ]
    SamAnnotator(annotate_file=args.tag_file,
                 samtags=samtags,
                 output_path=args.output_file,
                 allow_dovetailing=args.allow_dovetailing,
                 discard_bad_alt_tag=not args.keep_suboptimal_alternate_tags)

if __name__ == "__main__":
    main()