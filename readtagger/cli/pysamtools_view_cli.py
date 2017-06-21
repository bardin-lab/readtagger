import click
from readtagger.pysamtools_view import view
from readtagger import VERSION


@click.command()
@click.option('--input_bam',
              help='Alignment file from which to extract reads',
              required=True,
              type=click.Path(exists=True))
@click.option('--output_bam',
              help='Write extracted reads to this Alignment File',
              required=True,
              type=click.Path())
@click.option('--region',
              help='Extract reads from this region. Format is chr:start-end',
              required=True)
@click.version_option(version=VERSION)
def cli(**kwargs):
    """Extract reads from regin and write to new BAM alignment file."""
    return view(**kwargs)
