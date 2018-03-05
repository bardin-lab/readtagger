from readtagger.bwa import SimpleAligner


SEQ = 'ATGCATGCATCATGCATCGAGGGGATCATGCATGCATCATGCATCGAGGGGATCATGCATGCATCATGCATCGAGGGGATC'


def test_simplealigner(tmpdir):  # noqa: D103
    tmp = tmpdir.strpath
    with SimpleAligner(reference_sequences=SEQ, tmp_dir=tmp) as aligner:
        aligner.align(SEQ)
