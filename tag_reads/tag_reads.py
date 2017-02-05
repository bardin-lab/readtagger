import argparse
import pysam
import six


class SamTagProcessor(object):
    def __init__(self, source_path, tag_prefix_self, tag_prefix_mate, tag_mate=True):
        self.tag_mate = tag_mate
        self.tag_prefix_self = tag_prefix_self
        self.tag_prefix_mate = tag_prefix_mate
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
        tags = dict(ref=self.source_alignment.get_reference_name(r.tid),
                    pos=r.reference_start,
                    cigar=r.cigarstring,
                    sense='AS' if r.is_reverse else 'S')
        detail_tag_value = "R:{ref},POS:{pos},CIGAR:{cigar},S:{sense}".format(**tags)
        ref_tag_value = tags['ref']
        detail_tag = "%sD" % self.tag_prefix_self
        ref_tag = "%sR" % self.tag_prefix_self
        return [(detail_tag, detail_tag_value), (ref_tag, ref_tag_value)]

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
                    tag_d[r.query_name][fw_rev] = self.compute_tag(r)
                else:
                    tag_d[r.query_name] = {fw_rev: self.compute_tag(r)}
        return tag_d

    def add_mate(self):
        detail_tag = "%sD" % self.tag_prefix_mate
        ref_tag = "%sR" % self.tag_prefix_mate
        for qname, tag_d in six.iteritems(self.result):
            if len(tag_d) == 2:
                # Both mates aligned, add mate tag
                tag_d['r1'].extend([(detail_tag, tag_d['r2'][0][1]), (ref_tag, tag_d['r2'][1][1])])
                tag_d['r2'].extend([(detail_tag, tag_d['r1'][0][1]), (ref_tag, tag_d['r1'][1][1])])
            elif len(tag_d) == 1:
                # Only one of the mates mapped, so we fill the mate with its mate tag
                if 'r1' in tag_d:
                    tag_d['r2'] = [(detail_tag, tag_d['r1'][0][1]), (ref_tag, tag_d['r1'][1][1])]
                else:
                    tag_d['r1'] = [(detail_tag, tag_d['r2'][0][1]), (ref_tag, tag_d['r2'][1][1])]
            else:
                continue
            self.result[qname] = tag_d

    def get_tag(self, other_r):
        """
        convinience method that takes a read object `other_r` and fetches the
        annotation tag that has been processed in the instance
        """
        other_r_reverse = 'r1' if other_r.is_read1 else 'r2'
        tagged_mates = self.result.get(other_r.query_name)
        if tagged_mates:
            return tagged_mates[other_r_reverse]
        else:
            return None

class SamAnnotator(object):
    def __init__(self, annotate_file, samtags, output_path="test.bam"):
        """
        Compare `readtags` with `annotate_file`.
        Produces a new alignment file at output_path.
        :param annotate_file: 'Path to bam/sam file'
        :type annotate_file: str
        :param samtags: SamTagProcessor instance or list of SamTagProcessor instances
        :type samtags: SamTagProcessor or List[SamTagProcessor]
        """
        self.annotate_file = pysam.AlignmentFile(annotate_file)
        self.output = pysam.AlignmentFile(output_path, 'wb', template=self.annotate_file)
        self.samtags = samtags
        if isinstance(self.samtags, list):
            self.process_samtag_list()
        else:
            self.process()

    def process(self):
        """
        Walk along reads in self.annotate_file, add tags if read is stored in self.samtags
        and write out alignment
        """
        for read in self.annotate_file:
            annotated_tag = self.samtags.get_tag(read)
            if annotated_tag:
                read.tags = annotated_tag + read.tags
            self.output.write(read)
        self.output.close()

    def process_samtag_list(self):
        for read in self.annotate_file:
            for samtag in self.samtags:
                annotated_tag = samtag.get_tag(read)
                if annotated_tag:
                    read.tags = annotated_tag + read.tags
            self.output.write(read)
        self.output.close()


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
        if six.PY3:
            filetag = filetag.decode()  # py3 bytestring to unicode
        if ':' in filetag:
            filepath, tag, tag_mate = filetag.split(':')
            annotate_with.append(filepath)
            tag_prefix.append(tag.upper())
            tag_prefix_mate.append(tag_mate.upper())
        else:
            annotate_with.append(filetag)
            tag_prefix.append('R')  # Default is R for read, M for mate
            tag_prefix_mate.append('M')
    return annotate_with, tag_prefix, tag_prefix_mate



def parse_args():
    p = argparse.ArgumentParser(description="Tag reads in an alignment file based on other alignment files",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-t', '--tag_file', help="Tag reads in this file.", required=True)
    p.add_argument('-a', '--annotate_with', help="Tag reads in readfile if reads are aligned in these files."
                                                 "Append `:A:B` to tag first letter of tag describing read as A, "
                                                 "and first letter of tag describing the mate as B",
                   nargs = "+",
                   required=True)
    p.add_argument('-o', '--output_file', help="Write bam file to this path", required=True)
    return p.parse_args()

def main():
    args = parse_args()
    files_tags = zip(*parse_file_tags(args.annotate_with))
    samtags = [SamTagProcessor(filepath, tag_prefix_self=tag, tag_prefix_mate=tag_mate) for (filepath, tag, tag_mate) in files_tags ]
    SamAnnotator(annotate_file=args.tag_file, samtags=samtags, output_path=args.output_file)

if __name__ == "__main__":
    main()
