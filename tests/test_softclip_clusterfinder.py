from readtagger.find_softclip_clusters import (
    SoftClipClusterFinder,
)

WRONG_TSD = 'wrong_tsd.bam'
MULTI_H6 = 'multisample_h6.bam'


def test_softclip_clusterfinder(datadir_copy, tmpdir):  # noqa: D103, F811
    input_path = str(datadir_copy[WRONG_TSD])
    gff_out = tmpdir.join('out.gff').strpath
    clusterfinder = SoftClipClusterFinder(input_path=input_path, output_gff=gff_out, sample_name='sample')
    clusters = clusterfinder.clusters
    assert len(clusters) == 8


def test_softclip_clusterfinder_multi_h6(datadir_copy, tmpdir):  # noqa: D103, F811
    input_path = str(datadir_copy[MULTI_H6])
    gff_out = tmpdir.join('out.gff').strpath
    clusterfinder = SoftClipClusterFinder(input_path=input_path, output_gff=gff_out, sample_name='sample')
    clusters = clusterfinder.clusters
    assert len(clusters) == 8
