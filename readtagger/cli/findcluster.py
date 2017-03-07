import click
from readtagger import findcluster
from readtagger.readtagger import __VERSION__


@click.command()
@click.option('--input_path',
              help='Find cluster in this BAM file.',
              type=click.Path(exists=True))
@click.option('--output_bam',
              help='Write out BAM file with cluster information to this path. '
                   'Reads will have an additional "CD" tag to indicate the cluster number',
              type=click.Path(exists=False))
@click.option('--output_gff',
              help='Write out GFF file with cluster information to this path.',
              type=click.Path(exists=False))
@click.option('--sample_name',
              default=None,
              help='Sample name to use when writing out clusters in GFF file. '
                   'Default is to infer the name from the input filename.',
              )
@click.option('--include_duplicates/--no-include_duplicates',
              help='Include reads marked as duplicates when finding clusters.',
              default=False)
@click.option('--threads', help='Threads to use for cap3 assembly step', default=1, type=click.IntRange(1, 100))
@click.version_option(version=__VERSION__)
def cli(**kwds):
    """Find clusters of reads that support a TE insertion."""
    return findcluster.ClusterFinder(**kwds)
