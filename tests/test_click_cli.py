from readtagger.cli.allow_dovetailing import allow_dovetailing
from readtagger.cli.add_matesequence import annotate_mate
from readtagger.cli.annotate_softclipped_reads import annotate_softclipped_reads
from readtagger.cli.findcluster import findcluster
from readtagger.cli.plot_coverage import plot_coverage
from readtagger.cli.pysamtools_view_cli import pysamtools_view
from readtagger.cli.update_mapq import update_mapq
from readtagger.cli.write_supplementary_fastq import write_supplementary_fastq
from readtagger import VERSION
from click.testing import CliRunner
import pytest
from .helpers import reference_fasta  # noqa: F401

command_line_functions = [allow_dovetailing,
                          annotate_mate,
                          annotate_softclipped_reads,
                          findcluster,
                          plot_coverage,
                          pysamtools_view,
                          update_mapq,
                          write_supplementary_fastq]
ids = [f.name for f in command_line_functions]
EXTENDED = 'extended_annotated_updated_all_reads.bam'


@pytest.mark.parametrize('f', command_line_functions, ids=ids)
def test_click_cli_help(f):  # noqa: D103
    runner = CliRunner()
    result = runner.invoke(f, ['--help'])
    assert result.exit_code == 0
    assert 'Usage:' in result.output


@pytest.mark.parametrize('f', command_line_functions, ids=ids)
def test_click_cli_no_option(f):  # noqa: D103
    runner = CliRunner()
    result = runner.invoke(f)
    assert result.exit_code == 2
    assert 'Error: Missing' in result.output


@pytest.mark.parametrize('f', command_line_functions, ids=ids)
def test_click_cli_get_version(f):  # noqa: D103
    runner = CliRunner()
    result = runner.invoke(f, ['--version'])
    assert result.exit_code == 0
    assert VERSION in result.output


def test_add_matesequence_cli(datadir_copy, tmpdir):  # noqa: D103
    in_path = str(datadir_copy[EXTENDED])
    out_path = tmpdir.join('out.bam').strpath
    runner = CliRunner()
    result = runner.invoke(annotate_mate, ['--source', in_path, '--target', in_path, '--output_path', out_path])
    assert result.exit_code == 0


def test_annotate_softclipped_reads_cli(datadir_copy, tmpdir, reference_fasta):  # noqa: D103,F811
    in_path = str(datadir_copy[EXTENDED])
    out_path = tmpdir.join('out.bam').strpath
    runner = CliRunner()
    result = runner.invoke(annotate_softclipped_reads, ['--source',
                                                        in_path,
                                                        '--reference_fasta',
                                                        reference_fasta,
                                                        '--output_path',
                                                        out_path])
    assert result.exit_code == 0


def test_pysamtools_view_cli(datadir_copy, tmpdir):  # noqa: D103
    in_path = str(datadir_copy[EXTENDED])
    out_path = tmpdir.join('out.bam').strpath
    runner = CliRunner()
    result = runner.invoke(pysamtools_view, ['--input_bam',
                                             in_path,
                                             '--region',
                                             '2L:10-100',
                                             '--output_bam',
                                             out_path])
    assert result.exit_code == 0
