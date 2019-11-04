import click

from readtagger.summarize_fragments import (
    summarize_reads,
    write_evidence,
)


@click.command("summarize_reads")
@click.argument('input_path', type=click.Path(exists=True))
@click.argument('output_path', type=click.Path())
def cli(input_path, output_path):  # noqa: D103
    summary = summarize_reads(input_path)
    write_evidence(summary, output_path)
