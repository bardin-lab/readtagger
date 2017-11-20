import click

from readtagger.plot_coverage import plot_coverage_in_regions


@click.command()
@click.argument('file1')
@click.argument('label1')
@click.argument('file2')
@click.argument('label2')
@click.argument('output_path')
@click.option('--cores', default=1)
@click.option('--regions')
def script_entry(**kwargs):
    """Plot coverage differences between file1 and file2."""
    files = [kwargs.pop('file1'), kwargs.pop('file2')]
    labels = [kwargs.pop('label1'), kwargs.pop('label2')]
    regions = kwargs.get('regions')
    if regions:
        kwargs['regions'] = regions.split(',')
    plot_coverage_in_regions(files, labels, **kwargs)


if __name__ == '__main__':
    script_entry()
