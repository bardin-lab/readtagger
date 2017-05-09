import click

from readtagger import VERSION
from readtagger.tag_softclip import TagSoftClip


@click.command()
@click.option('--source',
              help='Alignment file for which to search additional alignments for softclipped reads',
              required=True,
              type=click.Path(exists=True))
@click.option('--reference_fasta',
              help='Fasta file to align softclipped reads against',
              required=True,
              type=click.Path(exists=True))
@click.option('--output_path',
              help='Write reads with updated tags here',
              required=True,
              type=click.Path(exists=False))
@click.option('--threads',
              help='Number of threads to use in BWA alignment',
              default=1)
@click.option('--min_clip_length',
              help='Minimum length for a clipped region to be extracted',
              default=20)
@click.version_option(version=VERSION)
def tag(**kwargs):
    """Dispatch CLI arguments to TagSoftClip."""
    TagSoftClip(**kwargs)
