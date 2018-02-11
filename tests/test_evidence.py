import pysam

from readtagger.cluster import evidence_for

TEST_BAM = 'improve_counting_support.bam'
GOOD_BREAKPOINTS = [13813889, 13813894]
WRONG_BREAKPOINTS = [13813890]
WRONG_BREAKPOINT_SEQS = ['ATATA', 'ATATA']
GOOD_BREAKPOINT_SEQS = ['TAAGT', 'GTAAC']


def test_evidence_for(datadir_copy):  # noqa: D103
    three_p_qname = 'HISEQ:788:H2G5VBCX2:1:2109:14202:61591'
    five_p_qname = 'HISEQ:788:H2G5VBCX2:1:2109:11125:90225'
    random_qname = 'HISEQ:785:H2G75BCX2:2:1202:6302:39767'

    three_p_read = _get_read(datadir_copy, qname=three_p_qname, is_read1=False)
    _breakpoint_and_sequence_combinations(three_p_read, should_be_true=True)

    five_p_read = _get_read(datadir_copy, qname=five_p_qname, is_read1=True)
    _breakpoint_and_sequence_combinations(five_p_read, should_be_true=True)

    random_read = _get_read(datadir_copy, qname=random_qname, is_read1=True)
    _breakpoint_and_sequence_combinations(random_read, should_be_true=False)


def _get_read(datadir_copy, qname, is_read1):
    p = str(datadir_copy[TEST_BAM])
    with pysam.AlignmentFile(p) as f:
        reads = [r for r in f if r.query_name == qname and r.is_read1 == is_read1]
        assert len(reads) == 1, "Should have 1 read, but got %d reads" % len(reads)
        return reads[0]


def _breakpoint_and_sequence_combinations(read, should_be_true=True):
    assert evidence_for(read=read, breakpoints=GOOD_BREAKPOINTS, breakpoint_sequences=WRONG_BREAKPOINT_SEQS) is False, str(read)
    assert evidence_for(read=read, breakpoints=WRONG_BREAKPOINTS, breakpoint_sequences=GOOD_BREAKPOINT_SEQS) is False, str(read)
    assert evidence_for(read=read, breakpoints=GOOD_BREAKPOINTS, breakpoint_sequences=GOOD_BREAKPOINT_SEQS) == should_be_true, str(read)
