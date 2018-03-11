import click

from readtagger.plot_coverage import plot_coverage_in_regions
from readtagger import VERSION


@click.command('Plot relative coverage for alignment files.')
@click.argument('file1')
@click.argument('label1')
@click.argument('file2')
@click.argument('label2')
@click.argument('output_path')
@click.option('-c',
              '--cores',
              help='Cores to use when calculating coverage',
              default=1)
@click.option('-r',
              '--regions',
              help='Regions to plot. If not specified plots all contigs.')
@click.version_option(version=VERSION)
def plot_coverage(**kwargs):
    """Plot coverage differences between file1 and file2."""
    files = [kwargs.pop('file1'), kwargs.pop('file2')]
    labels = [kwargs.pop('label1'), kwargs.pop('label2')]
    regions = kwargs.get('regions')
    if regions:
        kwargs['regions'] = regions.split(',')
    plot_coverage_in_regions(files, labels, **kwargs)
