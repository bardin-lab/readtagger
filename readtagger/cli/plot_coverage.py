import click

from readtagger.plot_coverage import plot_coverage_in_regions
from readtagger import VERSION


@click.command('Plot relative coverage for alignment files.')
@click.option('-f',
              '--file',
              type=(str, str, int),
              multiple=True,
              help="File, label and number of total reads in file.")
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
    file_tuples = kwargs.pop('file')
    kwargs['files'] = [_[0] for _ in file_tuples]
    kwargs['labels'] = [_[1] for _ in file_tuples]
    kwargs['total_reads'] = [_[2] for _ in file_tuples]
    regions = kwargs.get('regions')
    if regions:
        kwargs['regions'] = regions.split(',')
    plot_coverage_in_regions(**kwargs)
