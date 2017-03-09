from .bam_io import (
    BamAlignmentReader as Reader,
)


def write_supplementary_fastq(input_path, output_path):
    """
    Find supplementary reads in `input_path` and write them out to a fastq file at `output_path`.

    If supplementary reads are softclipped, only write out the matched sequence.
    """
    with Reader(input_path) as reader, open(output_path, 'w') as writer:
        for r in reader:
            if r.is_supplementary:
                seq = r.query_alignment_sequence
                qual = r.qqual
                name = r.query_name
                writer.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))
