from readtagger.findcluster import ClusterFinder

NON_SUPPORT = 'non_support_test.bam'


def test_nonevidence(datadir_copy):  # noqa: D103
    input_path = str(datadir_copy[NON_SUPPORT])
    result = ClusterFinder(input_path=input_path).clusters[:4]
    assert len(result) == 4
    assert result[0].nref == 12
    assert result[1].nref == 39
    assert result[2].nref == 30
    assert result[3].nref == 32
