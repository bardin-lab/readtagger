import click

from readtagger.extract_variants import extract_variant_portion


@click.command(name='Extract inserts and softclipped alignments')
@click.argument('input_path', nargs=1, type=click.Path(exists=True), required=True)
@click.argument('output_path', nargs=1, required=True)
@click.argument('filter_variant_source', nargs=1, required=False)
def extract_variants(input_path, output_path, filter_variant_source):
    """Extract unaligned part of alignment in ``input_path`` and write to ``output_path``."""
    extract_variant_portion(input_path, output_path, filter_variant_source)
