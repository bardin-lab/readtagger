import click
from readtagger.update_mapq import update_mapq
from readtagger import VERSION


@click.command()
@click.option('--source_path',
              help='Alignment file where supplementary reads should be updated with new MAPQ scores',
              required=True,
              type=click.Path(exists=True))
@click.option('--remapped_path',
              help='Alignment file that contains reads that have been remapped and contain an updated MAPQ score',
              required=True,
              type=click.Path(exists=True))
@click.option('--output_path',
              help='Write all alignments in original_path to this location. Supplementary reads will take the MAPQ score '
                   'as determined by remapping to the alignment file at `remapped_path`',
              required=True,
              type=click.Path(exists=False))
@click.version_option(version=VERSION)
def cli(**kwargs):
    """Update supplementary read MAPQ score after remapping."""
    return update_mapq(**kwargs)
