from click.testing import CliRunner
from readtagger.cli.allow_dovetailing import allow_dovetailing
from readtagger.allow_dovetailing import (
    get_max_proper_pair_size,
)
from readtagger.bam_io import BamAlignmentReader as Reader

NOT_A_PAIR = 'not_a_proper_pair.bam'
EXTENDED = 'extended_annotated_updated_all_reads.bam'


def test_allow_dovetailing(datadir_copy):  # noqa: D103
    max_proper_size = get_max_proper_pair_size(str(datadir_copy[EXTENDED]), reads_to_check=2)
    assert max_proper_size == 190
    max_proper_size = get_max_proper_pair_size(str(datadir_copy[EXTENDED]), reads_to_check=1000)
    assert max_proper_size == 664
    max_proper_size = get_max_proper_pair_size(str(datadir_copy[NOT_A_PAIR]), reads_to_check=1)
    assert not max_proper_size


def test_main(datadir_copy, tmpdir, mocker):  # noqa: D103
    outpath = tmpdir.join('pair.bam').strpath
    input_path = str(datadir_copy[NOT_A_PAIR])
    runner = CliRunner()
    result = runner.invoke(allow_dovetailing, ['--input_path', input_path, '--output_path', outpath])
    assert result.exit_code == 0
    with Reader(outpath) as reader:
        assert len([r for r in reader if r.is_proper_pair]) == 2
