import click
from readtagger.write_supplementary_fastq import write_supplementary_fastq as _write_supplementary_fastq
from readtagger import VERSION


@click.command()
@click.option('-i',
              '--input_path',
              help='Collect supplementary reads in this alignment file.',
              required=True,
              type=click.Path(exists=True))
@click.option('-o',
              '--output_path',
              help='Write supplementary reads to this FASTQ file.',
              required=True,
              type=click.Path(exists=False))
@click.version_option(version=VERSION)
def write_supplementary_fastq(input_path, output_path):
    """Write all supplementary alignments in `input_path` to an output fastq file at `output_path`."""
    return _write_supplementary_fastq(input_path=input_path, output_path=output_path)
