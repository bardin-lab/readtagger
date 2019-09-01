import click
from readtagger import VERSION
from readtagger.filter_insertions import confirm_insertions as _confirm_insertions


@click.command()
@click.option('-p',
              '--putative_insertions_path',
              help='Path to file containing putative somatic insertions',
              required=True,
              type=click.Path(exists=True))
@click.option('-t',
              '--all_treatments_path',
              help='Path to file containing all treatment insertions',
              required=True,
              type=click.Path(exists=True))
@click.option('-c',
              '--all_controls_paths',
              help='Path to file containing all control insertions',
              required=True,
              multiple=True,
              type=click.Path(exists=True))
@click.option('-o',
              '--output_path',
              help='Write annotated output to this Path.')
@click.option('--min_length',
              help="Minimum length necessary to match clip patterns."
                   "If too low will create false matches, if too high will not find legitimate matches",
              default=4,
              type=click.INT,
              )
@click.option('--output_discarded_records/--no_output_discarded_records',
              default=True,
              help="Discard an alternative flag if the current read is in a proper pair.")
@click.version_option(version=VERSION)
def confirm_insertions(**kwargs):
    """Confirm insertions by checking that control file does not contain the same clipping pattern."""
    _confirm_insertions(**kwargs)
