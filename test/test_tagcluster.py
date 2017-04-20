from readtagger.bam_io import BamAlignmentReader as Reader
from readtagger.tagcluster import TagCluster
from readtagger.tagcluster import TargetSiteDuplication

INPUT = 'tagged_dm6.bam'


def test_tsd(datadir):  # noqa: D103
    cluster = get_cluster(datadir)
    tsd = TargetSiteDuplication(cluster)
    assert len(tsd.unassigned_support) == 0
    assert tsd.is_valid
    tsd = TargetSiteDuplication(cluster, include_duplicates=True)
    assert len(tsd.unassigned_support) == 1
    assert tsd.is_valid
    [r.set_tag('AD', None) for r in cluster]
    tsd = TargetSiteDuplication(cluster, include_duplicates=True)
    assert not tsd.is_valid


def test_tagcluster_with_splits(datadir):  # noqa: D103
    cluster = get_cluster(datadir)
    tc = TagCluster(cluster)
    assert len(tc.tsd.three_p_support) == 4
    assert len(tc.tsd.five_p_support) == 1
    assert len(tc.tsd.unassigned_support) == 0
    assert tc.joint_insert.assembly.ncontigs == 0


def test_tagcluster_without_splits(datadir):  # noqa: D103
    cluster = get_cluster(datadir)
    [r.set_tag('AD', None) for r in cluster]
    tc = TagCluster(cluster)
    assert tc.joint_insert.assembly.ncontigs == 1


def get_cluster(datadir):
    """Get readcluster(= all reads) from input."""
    with Reader(datadir[INPUT]) as reader:
        return [r for r in reader]
