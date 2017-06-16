import copy
import os
import tempfile

import pysam
from six import string_types

from .bwa import Bwa
from .cigar import (
    cigar_tuple_to_cigar_length,
)
from .tags import Tag

SOFT_CLIP = 4  # This is the cigar operation for softclipped reads


class TagSoftClip(object):
    """Extract all softclipped sequences in a set of reads."""

    def __init__(self, source, reference_fasta=None, bwa_index=None, output_path=None, threads=1, min_clip_length=20):
        """
        Extract softclipped positions for reads in source and align against refrence_fasta.

        :param source: path to alignment file or list of reads
        :param reference_fasta: Reference to align clipped portion of reads to.
                                Only used if bwa_index is not given.
        :param bwa_index: Reference to align clipped portion of reads to.
                          Takes precedence over reference_fasta.
        :param output_path: path to output alignment file or None
        :param min_clip_length: minimum length of clip portion to output
        """
        self.source = source
        self.reference_fasta = reference_fasta
        self.bwa_index = bwa_index
        self.output_path = output_path
        self.threads = threads
        self.min_clip_length = min_clip_length
        fd, self.fastq_tmp = tempfile.mkstemp()
        os.close(fd)
        self.setup()
        self.write_clipped_portion()
        self.bwa = self.align_clipped_portion()
        self.aligned_clipped_reads = self.process_aligned_clipped_reads()
        self.annotate_clipped_reads()

    def setup(self):
        """Setup input and output files if these are paths."""
        if isinstance(self.source, string_types):
            self.source = pysam.AlignmentFile(self.source)
        if isinstance(self.output_path, string_types):
            self.writer = pysam.AlignmentFile(self.output_path, mode='wb', header=self.source.header)

    def write_clipped_portion(self):
        """Generate a list of reads to annotate with softclip alternative match."""
        with open(self.fastq_tmp, 'w') as fastq:
            for read in self.source:
                if not read.is_unmapped:
                    softclipped_portions = self._get_softclipped_portion(read)
                    for read, start, end in softclipped_portions:
                        fastq.write(self.clip_to_fastq(read, start, end))
        if isinstance(self.source, pysam.AlignmentFile):
            self.source.reset()

    def _get_softclipped_portion(self, read):
        softclipped_portions = []
        cigar_lengths = cigar_tuple_to_cigar_length(read.cigar)
        if len(cigar_lengths) > 1:
            for (start, end), operation in (cigar_lengths[0], cigar_lengths[-1]):
                if operation == SOFT_CLIP and end - start >= self.min_clip_length:
                    softclipped_portions.append((read, start, end))
        return softclipped_portions

    def align_clipped_portion(self):
        """Align fasta file against reference fasta."""
        return Bwa(input_path=self.fastq_tmp, reference_fasta=self.reference_fasta, bwa_index=self.bwa_index, threads=self.threads, describe_alignment=False)

    def process_aligned_clipped_reads(self):
        """Read newly aligned clipped reads into memory."""
        annotated_clipped = {}
        for read in self.bwa.bwa_run:
            if not read.is_unmapped:
                query_name = read.query_name
                if query_name not in annotated_clipped:
                    annotated_clipped[read.query_name] = [read]
                else:
                    annotated_clipped[read.query_name].append(read)
        return annotated_clipped

    def annotate_clipped_reads(self):
        """Iterate over source reads and fetch an alternative clipped alignment from self.aligned_clipped_reads."""
        for read in self.source:
            if not read.is_unmapped:
                softclipped_portions = self._get_softclipped_portion(read)
                for read, start, end in softclipped_portions:
                    potential_qname = "%s|%s|%s" % (read.query_name, start, end)
                    if potential_qname in self.aligned_clipped_reads:
                        for annotated_read in self.aligned_clipped_reads[potential_qname]:
                            tag = self.make_tag(start=start, end=end, original_read=read, annotated_clipped_read=annotated_read)
                            if read.has_tag('AD'):
                                tag = "%s,%s" % (read.get_tag('AD'), tag)
                            read.set_tag('AD', tag)
            if isinstance(self.output_path, string_types):
                self.writer.write(read)

    def make_tag(self, start, end, original_read, annotated_clipped_read):
        """Make a tag from the clipped alignment that can be added to the original read."""
        cigar = copy.copy(annotated_clipped_read.cigar)
        if start == 0:
            # The clipped portion started at 0, so in terms of the original read
            # we should have an alignment in the beginning of the read and a softclipped
            # portion at the end of the read (which is properly aligned in the original read).
            # This means that if reads have a softclipped portions that laigns to another reference,
            # the clipped parts should be complementary, i.e 82S43M and 43S82M.
            new_end = original_read.query_length - end
            cigar.append((SOFT_CLIP, new_end))
        else:
            cigar.insert(0, (SOFT_CLIP, start))
        if original_read.is_reverse != annotated_clipped_read.is_reverse:
            cigar = cigar[::-1]
        return Tag(tid=annotated_clipped_read.tid,
                   reference_start=annotated_clipped_read.reference_start,
                   cigar=cigar,
                   is_reverse=annotated_clipped_read.is_reverse,
                   mapq=annotated_clipped_read.mapping_quality,
                   query_alignment_start=annotated_clipped_read.query_alignment_start,
                   query_alignment_end=annotated_clipped_read.query_alignment_end,
                   ).to_string(header=self.bwa.header)

    def clip_to_fastq(self, read, start, end):
        """Select relevant portion of read to output as fastq."""
        seq_to_check = read.query_sequence[start:end]
        qual = self.phred_qual_to_sanger_string(read.query_qualities[start:end])
        return "@%s|%s|%s\n%s\n+\n%s\n" % (read.query_name, start, end, seq_to_check, qual)

    @staticmethod
    def phred_qual_to_sanger_string(qual):
        """Transform phred-scaled quality values to sanger string representation."""
        return "".join([chr(x + 33) for x in qual])
