import click

from readtagger.create_multisample_vcf import VCFMerger
from readtagger import VERSION


@click.command(name='Merge findcluster VCF files')
@click.argument('input_files', nargs=-1, type=click.Path(exists=True), required=True)
@click.argument('output_path', nargs=1, required=True)
@click.version_option(version=VERSION)
def merge_findcluster(**kwargs):
    """Merge many individual VCF files produced by findcluster."""
    VCFMerger(variant_file_paths=kwargs.pop('input_files'), output_path=kwargs.pop('output_path'))
