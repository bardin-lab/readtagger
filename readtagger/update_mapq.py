from .bam_io import BamAlignmentReader as Reader
from .bam_io import BamAlignmentWriter as Writer


def update_mapq(source_path, remapped_path, output_path):
    """
    Update supplementary read MAPQ score after remapping.

    :param source_path: Alignment file where supplementary reads should be updated with new MAPQ scores
    :param remapped_path: Alignment file that contains reads that have been remapped and contain an updated MAPQ score
    :param output_path: Write all alignments in original_path to this location. Supplementary reads will take the MAPQ score
                        as determined by remapping to the alignment file at `remapped_path`
    """
    with Reader(source_path) as source, Reader(remapped_path) as remapped, Writer(output_path, template=source) as output:
        qname_seq = {}
        for r in remapped:
            if r.query_name not in qname_seq:
                qname_seq[r.query_name] = {r.query_alignment_sequence: r.mapq}
            else:
                qname_seq[r.query_name][r.query_alignment_sequence] = r.mapq
        for r in source:
            if r.query_name in qname_seq and r.query_alignment_sequence in qname_seq[r.query_name]:
                r.mapq = qname_seq[r.query_name][r.query_alignment_sequence]
            output.write(r)
