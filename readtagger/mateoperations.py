import pysam
from cached_property import cached_property
from six import string_types


class AnnotateMateInformation(object):
    """Reads along a file that has a complete set of alignments and a file that should be annotated for mates of interest."""

    def __init__(self, target, source, output_path=None, mate_sequence_tag='MS'):
        """
        Add mate sequence into MS tag.

        :param file_to_annotate: path to alignment file
        :param annotate_source: path to alignment file
        :param output_path: path to output alignment file
        """
        self.source = source
        self.target = target
        self.output_path = output_path
        self.mate_sequence_tag = mate_sequence_tag
        self.setup()
        self.reads_to_annotate = self.get_reads_to_annotate()
        self.get_mates()
        self.write_annotated_reads()

    def setup(self):
        """Setup input and output files if these are paths."""
        if isinstance(self.source, string_types):
            self.source = pysam.AlignmentFile(self.source)
        if isinstance(self.target, string_types):
            self.target = pysam.AlignmentFile(self.target)
        if isinstance(self.output_path, string_types):
            self.writer = pysam.AlignmentFile(self.output_path, mode='wb', header=self.target.header)

    @cached_property
    def header(self):
        """Return header of file to annotate."""
        return self.target.header

    def get_reads_to_annotate(self):
        """Generate a list of mates to annotate."""
        reads = {}
        for read in self.target:
            reads["%s_%d" % (read.query_name, int(read.is_read1))] = None
        if isinstance(self.target, pysam.AlignmentFile):
            self.target.reset()
        return reads

    def get_mates(self):
        """Iterate over source file and annotate self.reads_to_annotate with the needed information."""
        for read in self.source:
            mate_id = "%s_%d" % (read.query_name, int(not read.is_read1))
            if mate_id in self.reads_to_annotate:
                self.reads_to_annotate[mate_id] = read.query_sequence

    def write_annotated_reads(self):
        """Add mate sequence to read in input file and write out."""
        for read in self.target:
            read_id = "%s_%d" % (read.query_name, int(read.is_read1))
            mate_seq = self.reads_to_annotate[read_id]
            read.set_tag(self.mate_sequence_tag, mate_seq)
            if self.output_path:
                self.writer.write(read)
