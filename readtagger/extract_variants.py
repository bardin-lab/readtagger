import copy

import mappy as mp

from readtagger.bam_io import (
    BamAlignmentReader,
    BamAlignmentWriter,
    sort_bam,
)
from readtagger.cigar import (
    cigar_tuple_to_cigar_length,
    INSERTION,
    SOFT_CLIP,
)
from readtagger.tags import Tag


def extract_variant_portion(input_path, output_path, filter_variant_source=None):
    """
    Extract unaligned part of alignment in ``input_path`` and write to ``output_path``.

    The unaligned part may  be clipped sequences or insertions in the alignment.
    """
    if filter_variant_source:
        aligner = mp.Aligner(filter_variant_source, preset='map-ont')
    with BamAlignmentReader(input_path) as bam_in, BamAlignmentWriter(output_path, template=bam_in) as bam_out:
        for alignment in bam_in:
            for new_alignment_start, variant_sequence, cigar in yield_variant_sequences(alignment):
                if filter_variant_source:
                    hits = []
                    three_p_clip = variant_sequence.startswith('N')
                    for hit in aligner.map(variant_sequence.replace('N', '')):
                        if (three_p_clip and hit.q_st < 100) or (not three_p_clip and abs(hit.q_en - len(variant_sequence)) < 100):
                            hits.append(Tag(
                                reference_start=hit.r_st,
                                cigar=hit.cigar_str,
                                is_reverse=hit.strand == -1,
                                mapq=hit.mapq,
                                query_alignment_start=hit.q_st,
                                query_alignment_end=hit.q_en,
                                reference_name=hit.ctg
                            ).to_string())
                    if hits:
                        ad_tag = ";".join(hits)
                        # just pick the first contig for simplicity
                        ar_tag = hits[0].split(',')[0][2:]
                        tags = [('AD', ad_tag), ('AR', ar_tag)]
                        new_alignment = copy.deepcopy(alignment)
                        new_alignment.tags = tags
                        new_alignment.reference_start = new_alignment_start
                        new_alignment.cigar = cigar
                        new_alignment.query_sequence = variant_sequence
                        bam_out.write(new_alignment)
    sort_bam(inpath=output_path, output=output_path, sort_order='coordinate', threads=1)


def yield_variant_sequences(alignment):
    """Yield soft-clipped and inserted sequences."""
    for (start, end), operation in cigar_tuple_to_cigar_length(alignment.cigar):
        if end - start > 100 and alignment.query_sequence and operation in (INSERTION, SOFT_CLIP):
            variant_sequence = alignment.query_sequence[start:end]
            if start == 0:
                # 5 prime soft clipping
                if end >= 100:
                    yield alignment.reference_start, "%sN" % variant_sequence, [(SOFT_CLIP, len(variant_sequence)), (0, 1)]
            elif end == alignment.query_length:
                yield alignment.reference_end, "N%s" % variant_sequence, [(0, 1), (SOFT_CLIP, len(variant_sequence))]
            else:
                # This has to be an insertion, need to infer start/end from aligned pairs
                ap = alignment.get_aligned_pairs()
                last_reference_position = None
                reference_start = None
                for (query_position, reference_position) in ap:
                    if reference_start is None and query_position and (query_position == start - 1 or query_position >= start):
                        # This should be the last aligned position, since deletions are None
                        reference_start = reference_position or last_reference_position
                        yield reference_start, "N%sN" % variant_sequence, [(0, 1), (INSERTION, len(variant_sequence)), (0, 1)]
                        break
                    if reference_position is not None:
                        last_reference_position = reference_position
