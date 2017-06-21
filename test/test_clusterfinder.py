from collections import namedtuple
from readtagger.findcluster import (
    ClusterFinder,
    ClusterManager
)
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
GENOME_FRAGMENT = 'genome_fragment.fa'
COMPLEX = 'extended_annotated_updated_all_reads.bam'
REORGANIZE_CLUSTER = 'reorganize_cluster.bam'
NON_SUPPORT = 'non_support_test.bam'
REFINE_COORD = 'refine_coord.bam'
REASSEMBLE = 'reassemble_input.bam'
SPLIT_CLUSTER = 'hum3_false_merge.bam'
SPLIT_CLUSTER_OPT = 'split_cluster_opt.bam'
SPLIT_CLUSTER_OPT2 = 'improve_clustering.bam'
MULTIPROCESSING = 'pasteurianus.bam'


def test_clusterfinder_single_cluster(datadir):  # noqa: D103
    input_path = datadir[INPUT]
    cf = ClusterFinder(input_path=input_path)
    assert len(cf.cluster) == 1
    assert cf.cluster[0].nalt == 19


def test_clusterfinder_include_duplicates(datadir):  # noqa: D103
    input_path = datadir[INPUT]
    cf = ClusterFinder(input_path=input_path, include_duplicates=True)
    assert len(cf.cluster) == 1
    assert cf.cluster[0].nalt == 26


def test_clusterfinder_remove_supplementary(datadir):  # noqa: D103
    input_path = datadir[INPUT]
    cf = ClusterFinder(input_path=input_path, include_duplicates=True, remove_supplementary_without_primary=True)
    assert len(cf.cluster) == 1
    assert cf.cluster[0].nalt == 24


def test_clusterfinder_reassemble(datadir, tmpdir, reference_fasta):   # noqa: D103, F811
    input_path = datadir[REASSEMBLE]
    genome_reference_fasta = datadir[GENOME_FRAGMENT]
    output_bam = tmpdir.join('output.bam').strpath
    clusters = ClusterFinder(input_path=input_path,
                             output_bam=output_bam,
                             output_gff=tmpdir.join('output.gff').strpath,
                             genome_reference_fasta=genome_reference_fasta,
                             transposon_reference_fasta=reference_fasta,
                             max_proper_pair_size=480)
    assert len(clusters.cluster) == 1
    assert clusters.cluster[0].valid_tsd


def test_clusterfinder_split_cluster(datadir, tmpdir):  # noqa: D103
    input_path = datadir[SPLIT_CLUSTER]
    cf = ClusterFinder(input_path=input_path, output_bam=tmpdir.join('output.bam').strpath,
                       include_duplicates=False,
                       remove_supplementary_without_primary=False,
                       output_gff=tmpdir.join('out.gff').strpath)
    assert len(cf.cluster) == 3
    assert len(cf.cluster[1]) == 39
    assert len(cf.cluster[2]) == 27


def test_clusterfinder_refine_split(datadir, tmpdir):  # noqa: D103
    input_path = datadir[SPLIT_CLUSTER_OPT]
    output_bam = tmpdir.join('output.bam').strpath
    clusters = ClusterFinder(input_path=input_path,
                             output_bam=output_bam,
                             output_gff=tmpdir.join('output.gff').strpath,
                             transposon_reference_fasta=None,
                             max_proper_pair_size=480)
    cluster = clusters.cluster[0]
    genotype = cluster.genotype_likelihood()
    assert genotype.nref == 72
    assert genotype.nalt == 107
    assert genotype.genotype == 'heterozygous'


def test_clusterfinder_refine_split2(datadir, tmpdir):  # noqa: D103
    input_path = datadir[SPLIT_CLUSTER_OPT2]
    output_bam = tmpdir.join('output.bam').strpath
    clusters = ClusterFinder(input_path=input_path,
                             output_bam=output_bam,
                             output_gff=tmpdir.join('output.gff').strpath,
                             transposon_reference_fasta=None,
                             max_proper_pair_size=480)
    assert len(clusters.cluster) == 2
    cluster_one, cluster_two = clusters.cluster
    assert cluster_one.nalt == 2
    assert cluster_two.nalt == 134


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
    assert len(cf.cluster) == 2


def test_clusterfinder_multiple_cluster_gff(datadir, tmpdir):  # noqa: D103
    input_path = datadir[EXTENDED]
    output_gff = tmpdir.join('output.gff')
    cf = ClusterFinder(input_path=input_path, output_gff=output_gff.strpath)
    assert len(cf.cluster) == 2


def test_clusterfinder_cache_threads(datadir, tmpdir, mocker):  # noqa: D103
    input_path = datadir[EXTENDED]
    output_gff = tmpdir.join('output.gff')
    cf = ClusterFinder(input_path=input_path, output_gff=output_gff.strpath, threads=2)
    assert len(cf.cluster) == 2
    mocker.spy(cf.cluster[0], 'can_join')
    mocker.spy(cf.cluster[0], '_can_join')
    cf.join_clusters()
    assert cf.cluster[0].can_join.call_count == 4
    assert cf.cluster[0]._can_join.call_count == 0
    assert len(cf.cluster) == 2


def test_clustermanager_single_core(datadir, tmpdir):  # noqa: D103
    input_path = datadir[EXTENDED]
    output_gff = tmpdir.join('output.gff').strpath
    output_bam = tmpdir.join('output.bam').strpath
    ClusterManager(input_path=input_path, genome_reference_fasta=None, transposon_reference_fasta=None, output_bam=output_bam, output_gff=output_gff, threads=1)


def test_clustermanager_multiprocessing(datadir, tmpdir):  # noqa: D103
    input_path = datadir[MULTIPROCESSING]
    output_gff = tmpdir.join('output.gff').strpath
    output_bam = tmpdir.join('output.bam').strpath
    ClusterManager(input_path=input_path, genome_reference_fasta=None, transposon_reference_fasta=None, output_bam=output_bam, output_gff=output_gff, threads=2)


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
    args_template = namedtuple('ArgumentParser', 'input_path output_gff output_bam transposon_reference_fasta')
    args = args_template(input_path=input_path, output_bam=output_bam, output_gff=output_gff, transposon_reference_fasta=reference_fasta)
    argv = namedtuple_to_argv(args)
    mocker.patch('sys.argv', argv)
    mocker.patch('sys.exit')
    findcluster.cli()


def test_clusterfinder_complex_genotype(datadir, tmpdir, reference_fasta):  # noqa: D103, F811
    input_path = datadir[COMPLEX]
    output_bam = tmpdir.join('output.bam').strpath
    output_gff = tmpdir.join('output.gff').strpath
    output_fasta = tmpdir.join('output.fasta').strpath
    clusters = ClusterFinder(input_path=input_path,
                             output_bam=output_bam,
                             output_gff=output_gff,
                             transposon_reference_fasta=reference_fasta,
                             output_fasta=output_fasta)
    cluster = clusters.cluster[0]
    assert cluster.nalt == 20
    genotype = cluster.genotype_likelihood()
    assert genotype.nref == 17
    assert genotype.nalt == 20
    assert genotype.genotype == 'heterozygous'
    assert len(open(output_fasta).readlines()) == 6


def test_clusterfinder_nonsupport(datadir, tmpdir):  # noqa: D103
    input_path = datadir[NON_SUPPORT]
    output_bam = tmpdir.join('output.bam').strpath
    clusters = ClusterFinder(input_path=input_path, output_bam=output_bam, output_gff=None, transposon_reference_fasta=None)
    cluster = clusters.cluster[-1]
    assert cluster.nref == 25  # Could also be 26 -- need to figure that out.
    genotype = cluster.genotype_likelihood()
    assert genotype.nref == 25
    assert genotype.nalt == 1
    assert genotype.genotype == 'reference'


def test_clusterfinder_refine_coord(datadir, tmpdir):  # noqa: D103
    input_path = datadir[REFINE_COORD]
    output_bam = tmpdir.join('output.bam').strpath
    clusters = ClusterFinder(input_path=input_path,
                             output_bam=output_bam,
                             output_gff=tmpdir.join('output.gff').strpath,
                             transposon_reference_fasta=None,
                             max_proper_pair_size=480)
    cluster = clusters.cluster[-1]
    genotype = cluster.genotype_likelihood()
    assert genotype.nref == 0
    assert genotype.nalt == 5
    assert genotype.genotype == 'homozygous'


def test_clusterfinder_reorganize_cluster(datadir, tmpdir, reference_fasta):  # noqa: D103, F811
    input_path = datadir[REORGANIZE_CLUSTER]
    output_gff = tmpdir.join('output.gff').strpath
    clusters = ClusterFinder(input_path=input_path, output_bam=None, output_gff=output_gff, transposon_reference_fasta=reference_fasta)
    cluster = clusters.cluster[-1]
    genotype = cluster.genotype_likelihood()
    assert genotype.genotype == 'homozygous'
