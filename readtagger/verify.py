"""Verify that supplementary alignments without primary alignments are discarded."""
from .bam_io import (
    BamAlignmentReader as Reader,
    BamAlignmentWriter as Writer
)


def discard_supplementary(input_path, output_path):
    """
    Discard all supplementary alignments without primary alignments.

    Iterate over file at `input_path` and record all primary reads.
    Iterate again over `input_path`  and write all reads to `output_path`,
    except for reads that do not have a primary alignment in `input_path`
    """
    with Reader(input_path) as input_alignemnt:
        primary_alignments = set()
        [primary_alignments.add(r.query_name) for r in input_alignemnt if not r.is_supplementary]
    with Reader(input_path) as reader, Writer(output_path, template=reader) as writer:
        [writer.write(r) for r in reader if not r.is_supplementary or r.is_supplementary and r.query_name in primary_alignments]
