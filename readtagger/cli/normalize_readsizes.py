import click

from readtagger import VERSION
from readtagger.normalization import split_fastq_files


@click.command()
@click.option('-i', '--input_paths', required=True, type=click.Path(), multiple=True)
@click.option('-o', '--output_paths', required=True, type=click.Path(), multiple=True)
@click.version_option(version=VERSION)
def normalize_readsizes(input_paths, output_paths):
    """Normalize length of multiple long read fastq files."""
    split_fastq_files(input_paths, output_paths)
