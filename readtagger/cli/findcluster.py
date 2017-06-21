import click
from readtagger import findcluster
from readtagger import VERSION


@click.command()
@click.option('--input_path',
              help='Find cluster in this BAM file.',
              type=click.Path(exists=True))
@click.option('--region',
              help='Find clusters in this Region (Format is chrX:2000-1000).',
              default=None,)
@click.option('--max_proper_pair_size',
              help='Maximum proper pairs size. If not given will be inferred from the data.',
              default=0,)
@click.option('--output_bam',
              help='Write out BAM file with cluster information to this path. '
                   'Reads will have an additional "CD" tag to indicate the cluster number',
              type=click.Path(exists=False))
@click.option('--output_gff',
              help='Write out GFF file with cluster information to this path.',
              type=click.Path(exists=False))
@click.option('--output_fasta',
              help='Write out supporting evidence for clusters to this path.',
              type=click.Path(exists=False))
@click.option('--sample_name',
              default=None,
              help='Sample name to use when writing out clusters in GFF file. '
                   'Default is to infer the name from the input filename.',
              )
@click.option('--include_duplicates/--no-include_duplicates',
              help='Include reads marked as duplicates when finding clusters.',
              default=False)
@click.option('--transposon_reference_fasta',
              help=('Transposon fasta to align clipped reads to. '
                    'Not necessary if BWA index is provided.'),
              default=None,
              required=False)
@click.option('--transposon_bwa_index',
              help='Transposon BWA index to align clipped reads to',
              default=None,
              required=False)
@click.option('--genome_reference_fasta',
              help=('Genome fasta to align clipped reads to. '
                    'Not necessary if BWA index is provided.'),
              default=None,
              required=False)
@click.option('--genome_bwa_index',
              help='Genome BWA index to align clipped reads to',
              default=None,
              required=False)
@click.option('--threads', help='Threads to use for cap3 assembly step', default=1, type=click.IntRange(1, 100))
@click.option('--shm_dir', envvar="SHM_DIR", help='Path to shared memory folder', default=None, type=click.Path(exists=True))
@click.version_option(version=VERSION)
def cli(**kwds):
    """Find clusters of reads that support a TE insertion."""
    return findcluster.ClusterManager(**kwds)
