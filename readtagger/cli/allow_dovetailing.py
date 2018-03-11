import click
from readtagger.allow_dovetailing import process
from readtagger import VERSION


@click.command('allow_dovetailing')
@click.option('-i', '--input_path',
              help='Input alignment file to manipulate',
              required=True,
              type=click.Path(exists=True))
@click.option('-o', '--output_path',
              help='Output alignment file',
              required=True,
              type=click.Path(exists=False))
@click.version_option(version=VERSION)
def allow_dovetailing(**kwargs):
    """Dispatch to allow_dovetailing."""
    process(**kwargs)
