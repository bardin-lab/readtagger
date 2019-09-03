import pysam
from readtagger.cigar import aligned_segment_corresponds_to_transposable_element

TEST_DATA = 'aligned_segment_corresponds_to_transposable_element.bam'


def test_aligned_segment_corresponds_to_transposable_element(datadir_copy):  # noqa: D103, F811
    f = str(datadir_copy[TEST_DATA])
    with pysam.AlignmentFile(f) as af:
        r = next(iter(af))
    assert aligned_segment_corresponds_to_transposable_element(r)
