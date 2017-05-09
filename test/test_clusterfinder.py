from collections import namedtuple
from readtagger.findcluster import ClusterFinder
from readtagger.cli import findcluster

from .helpers import (  # noqa: F401
    namedtuple_to_argv,
    reference_fasta
)

INPUT = 'tagged_dm6.bam'
CORNERCASE = 'cornercase.bam'
CORNERCASE2 = 'cornercase2.bam'
CORNERCASE3 = 'cornercase3.bam'
EXTENDED = 'extended_and_annotated_roi.bam'
COMPLEX = 'extended_annotated_updated_all_reads.bam'
REORGANIZE_CLUSTER = 'reorganize_cluster.bam'
NON_SUPPORT = 'non_support_test.bam'


def test_clusterfinder_single_cluster(datadir):  # noqa: D103
    input_path = datadir[INPUT]
    cf = ClusterFinder(input_path=input_path)
    assert len(cf.cluster) == 1
    assert len(cf.cluster[0]) == 20


def test_clusterfinder_include_duplicates(datadir):  # noqa: D103
    input_path = datadir[INPUT]
    cf = ClusterFinder(input_path=input_path, include_duplicates=True)
    assert len(cf.cluster) == 1
    assert len(cf.cluster[0]) == 27


def test_clusterfinder_remove_supplementary(datadir):  # noqa: D103
    input_path = datadir[INPUT]
    cf = ClusterFinder(input_path=input_path, include_duplicates=True, remove_supplementary_without_primary=True)
    assert len(cf.cluster) == 1
    assert len(cf.cluster[0]) == 25


def test_cornercase(datadir, tmpdir):  # noqa: D103
    input_path = datadir[CORNERCASE]
    output_gff = tmpdir.join('output.gff').strpath
    cf = ClusterFinder(input_path=input_path, output_gff=output_gff, min_mapq=-1)
    assert len(cf.cluster) == 1
    input_path = datadir[CORNERCASE2]
    cf = ClusterFinder(input_path=input_path, output_gff=output_gff)
    assert len(cf.cluster) == 0
    input_path = datadir[CORNERCASE3]
    cf = ClusterFinder(input_path=input_path, output_gff=output_gff)


def test_clusterfinder_multiple_cluster(datadir, tmpdir):  # noqa: D103
    input_path = datadir[EXTENDED]
    output_bam = tmpdir.join('tagged_clusters.bam')
    cf = ClusterFinder(input_path=input_path, output_bam=output_bam.strpath)
    assert len(cf.cluster) == 3


def test_clusterfinder_multiple_cluster_gff(datadir, tmpdir):  # noqa: D103
    input_path = datadir[EXTENDED]
    output_gff = tmpdir.join('output.gff')
    cf = ClusterFinder(input_path=input_path, output_gff=output_gff.strpath)
    assert len(cf.cluster) == 3


def test_clusterfinder_cache_threads(datadir, tmpdir, mocker):  # noqa: D103
    input_path = datadir[EXTENDED]
    output_gff = tmpdir.join('output.gff')
    cf = ClusterFinder(input_path=input_path, output_gff=output_gff.strpath, threads=2)
    assert len(cf.cluster) == 3
    mocker.spy(cf.cluster[0], 'can_join')
    mocker.spy(cf.cluster[0], '_can_join')
    cf.join_clusters()
    assert cf.cluster[0].can_join.call_count == 4
    assert cf.cluster[0]._can_join.call_count == 0
    assert len(cf.cluster) == 3


def test_clusterfinder_multiple_cluster_gff_cli(datadir, tmpdir, mocker):  # noqa: D103
    input_path = datadir[EXTENDED]
    output_bam = tmpdir.join('output.bam').strpath
    output_gff = tmpdir.join('output.gff').strpath
    args_template = namedtuple('ArgumentParser', 'input_path output_gff output_bam')
    args = args_template(input_path=input_path, output_bam=output_bam, output_gff=output_gff)
    argv = namedtuple_to_argv(args)
    mocker.patch('sys.argv', argv)
    mocker.patch('sys.exit')
    findcluster.cli()


def test_clusterfinder_blast(datadir, tmpdir, mocker, reference_fasta):  # noqa: D103, F811
    input_path = datadir[EXTENDED]
    output_bam = tmpdir.join('output.bam').strpath
    output_gff = tmpdir.join('output.gff').strpath
    args_template = namedtuple('ArgumentParser', 'input_path output_gff output_bam reference_fasta')
    args = args_template(input_path=input_path, output_bam=output_bam, output_gff=output_gff, reference_fasta=reference_fasta)
    argv = namedtuple_to_argv(args)
    mocker.patch('sys.argv', argv)
    mocker.patch('sys.exit')
    findcluster.cli()


def test_clusterfinder_complex_genotype(datadir, tmpdir, reference_fasta):  # noqa: D103, F811
    input_path = datadir[COMPLEX]
    output_bam = tmpdir.join('output.bam').strpath
    output_gff = tmpdir.join('output.gff').strpath
    output_fasta = tmpdir.join('output.fasta').strpath
    clusters = ClusterFinder(input_path=input_path, output_bam=output_bam, output_gff=output_gff, reference_fasta=reference_fasta, output_fasta=output_fasta)
    cluster = clusters.cluster[0]
    assert len(cluster) == 20
    genotype = cluster.genotype_likelihood()
    assert genotype.nref == 17
    assert genotype.nalt == 20
    assert genotype.genotype == 'heterozygous'
    assert len(open(output_fasta).readlines()) == 6


def test_clusterfinder_nonsupport(datadir, tmpdir):  # noqa: D103
    input_path = datadir[NON_SUPPORT]
    output_bam = tmpdir.join('output.bam').strpath
    clusters = ClusterFinder(input_path=input_path, output_bam=output_bam, output_gff=None, reference_fasta=None)
    cluster = clusters.cluster[-1]
    assert cluster.nref == 25  # Could also be 26 -- need to figure that out.
    genotype = cluster.genotype_likelihood()
    assert genotype.nref == 25
    assert genotype.nalt == 1
    assert genotype.genotype == 'reference'


def test_clusterfinder_reorganize_cluster(datadir, tmpdir, reference_fasta):  # noqa: D103, F811
    input_path = datadir[REORGANIZE_CLUSTER]
    output_gff = tmpdir.join('output.gff').strpath
    clusters = ClusterFinder(input_path=input_path, output_bam=None, output_gff=output_gff, reference_fasta=reference_fasta)
    cluster = clusters.cluster[-1]
    genotype = cluster.genotype_likelihood()
    assert genotype.genotype == 'homozygous'
