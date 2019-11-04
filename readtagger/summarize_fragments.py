import collections

import pandas as pd
import pysam

from readtagger.tags import Tag

COLUMNS = [
    'reference',
    'start',
    'end',
    'insert_reference',
    'insert_start',
    'insert_end',
    'insert_orientation',
    'query_start',
    'query_end',
    'query_length',
    'read_length',
    'query_name',
    'evidence_type',
]


def summarize_reads(input_path):
    """Summarizes insert evidence for alignemnt file at ``input_path``."""
    evidence = collections.namedtuple('Evidence', " ".join(COLUMNS))
    evidence_by_reference = []
    with pysam.AlignmentFile(input_path) as af:
        for r in af:
            reference = r.get_tag('AR')
            ad = Tag.from_tag_str(r.get_tag('AD').split(';')[0])
            reference_start = ad.reference_start
            reference_end = reference_start + sum(t.length for t in ad.cigar if t.operation in (0, 2))
            qs = r.query_sequence
            if qs.startswith('N') and qs.endswith('N'):
                rtype = "insert"
            elif qs.startswith('N'):
                rtype = "3p_clip"
            else:
                rtype = "5p_clip"
            evidence_by_reference.append(evidence(
                r.reference_name,
                r.reference_start,
                r.reference_end,
                reference,
                reference_start,
                reference_end,
                "-" if ad.is_reverse else "+",
                ad.query_alignment_start,
                ad.query_alignment_end,
                ad.query_alignment_end - ad.query_alignment_start,
                r.query_length,
                r.query_name,
                rtype,
            ))
    evidence_by_reference = pd.DataFrame.from_records(evidence_by_reference)
    evidence_by_reference.columns = COLUMNS
    return evidence_by_reference


def write_evidence(summary, output_path):
    """Write pandas object to tabular file."""
    summary.to_csv(output_path, index=None, sep="\t")
