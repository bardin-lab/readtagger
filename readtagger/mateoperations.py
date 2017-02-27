import argparse

from cached_property import cached_property

from .bam_io import BamAlignmentReader as Reader
from .bam_io import BamAlignmentWriter as Writer
from .readtagger import __VERSION__


class AnnotateMateInformation(object):
    """Reads along a file that has a complete set of alignments and a file that should be annotated for mates of interest."""

    def __init__(self, file_to_annotate, annotate_source, output_path, mate_sequence_tag='MS'):
        """
        Add mate sequence into MS tag.

        :param file_to_annotate: path to alignment file
        :param annotate_source: path to alignment file
        :param output_path: path to output alignment file
        """
        self.file_to_annotate = file_to_annotate
        self.annotate_source = annotate_source
        self.output_path = output_path
        self.mate_sequence_tag = mate_sequence_tag
        self.reads_to_annotate = self.get_reads_to_annotate()
        self.get_mates()
        self.write_annotated_reads()

    @cached_property
    def header(self):
        """Return header of file to annotate."""
        with Reader(self.file_to_annotate, external_bin=None) as reader:
            return reader.header

    def get_reads_to_annotate(self):
        """Generate a list of mates to annotate."""
        reads = {}
        with Reader(self.file_to_annotate) as reader:
            for read in reader:
                reads["%s_%d" % (read.query_name, int(read.is_read1))] = None
        return reads

    def get_mates(self):
        """Iterate over source file and annotate self.reads_to_annotate with the needed information."""
        with Reader(self.annotate_source) as reader:
            for read in reader:
                mate_id = "%s_%d" % (read.query_name, int(not read.is_read1))
                if mate_id in self.reads_to_annotate:
                    self.reads_to_annotate[mate_id] = read.query_sequence

    def write_annotated_reads(self):
        """Add mate sequence to read in input file and write out."""
        with Reader(self.file_to_annotate) as reader, Writer(self.output_path, header=self.header) as writer:
            for read in reader:
                read_id = "%s_%d" % (read.query_name, int(read.is_read1))
                mate_seq = self.reads_to_annotate[read_id]
                read.set_tag(self.mate_sequence_tag, mate_seq)
                writer.write(read)


def parse_args():
    """Parse commandline arguments."""
    p = argparse.ArgumentParser(description='Annotate reads with Mate Sequence in tag field',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-a', '--file_to_annotate', help='Annotate reads in this file with their mate sequence', required=True)
    p.add_argument('-s', '--annotate_source', help='Use this file to find the mate sequence (can be same file as file_to_annotate)', required=True)
    p.add_argument('-o', '--output_path', help='Write resulting BAM file to this path', required=True)
    p.add_argument('-ms', '--mate_sequence_tag', help='Use this tag to store the mate sequence', default='MS')
    p.add_argument('--version', action='version', version=__VERSION__)
    return p.parse_args()


def main(args=None):
    """Annotate reads with Mate Sequence in tag field."""
    if not args:
        args = parse_args()
    AnnotateMateInformation(file_to_annotate=args.file_to_annotate,
                            annotate_source=args.annotate_source,
                            output_path=args.output_path,
                            mate_sequence_tag=args.mate_sequence_tag)


if __name__ == "__main__":
    main()
