from collections import namedtuple
from readtagger.cli.write_supplementary_fastq import cli
from .helpers import namedtuple_to_argv

TEST_BAM = 'supplementary.bam'
TEMPLATE_ARGS = namedtuple('args', ['input_path', 'output_path'])


def test_write_supplementary(datadir, tmpdir, mocker):  # noqa: D103
    out = tmpdir.join('out.fastq').strpath
    args = TEMPLATE_ARGS(input_path=datadir[TEST_BAM], output_path=out)
    argv = namedtuple_to_argv(args, 'write_supplementary_fastq.py')
    mocker.patch('sys.argv', argv)
    mocker.patch('sys.exit')
    cli()
    assert len(open(out).readlines()) == 36
