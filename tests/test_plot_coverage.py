from readtagger.plot_coverage import plot_coverage_in_regions
from readtagger.cli.plot_coverage import plot_coverage

BAM_FILE = 'split_cluster_opt.bam'


def test_plot_coverage(datadir_copy, tmpdir):  # noqa: D103
    bam = str(datadir_copy[BAM_FILE])
    pdf = tmpdir.join('plot_test.pdf').strpath
    plot_coverage_in_regions(files=[bam, bam], labels=['a', 'b'], output_path=pdf, regions=['3L'])


def test_plot_coverage_region_multicore(datadir_copy, tmpdir):  # noqa: D103
    bam = str(datadir_copy[BAM_FILE])
    pdf = tmpdir.join('plot_test.pdf').strpath
    plot_coverage_in_regions(files=[bam, bam], labels=['a', 'b'], regions=["3L:3280261-3281790"], output_path=pdf, cores=2)


def test_main(datadir_copy, tmpdir, mocker):  # noqa: D103
    outpath = tmpdir.join('test_plot_cli.pdf').strpath
    bam = str(datadir_copy[BAM_FILE])
    args = ['path_to_script', bam, 'bam1', bam, 'bam2', outpath, '--regions', '3L:3280261-3281890']
    mocker.patch('sys.argv', args)
    mocker.patch('sys.exit', lambda x: True)
    plot_coverage()
