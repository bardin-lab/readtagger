from readtagger.allow_dovetailing import (
    main,
    get_max_proper_pair_size,
)
from readtagger.bam_io import BamAlignmentReader as Reader
from collections import namedtuple
from .helpers import namedtuple_to_argv

NOT_A_PAIR = 'not_a_proper_pair.bam'
EXTENDED = 'extended_annotated_updated_all_reads.bam'


def test_allow_dovetailing(datadir):  # noqa: D103
    with Reader(datadir[EXTENDED], external_bin=None) as reader:
        max_proper_size = get_max_proper_pair_size(reader, reads_to_check=2)
        assert max_proper_size == 190
    with Reader(datadir[EXTENDED], external_bin=None) as reader:
        max_proper_size = get_max_proper_pair_size(reader, reads_to_check=1000)
        assert max_proper_size == 664
    with Reader(datadir[NOT_A_PAIR], external_bin=None) as reader:
        max_proper_size = get_max_proper_pair_size(reader, reads_to_check=1)
        assert not max_proper_size


def test_main(datadir, tmpdir, mocker):  # noqa: D103
    outpath = tmpdir.join('pair.bam')
    args_template = namedtuple('ArgumentParser', 'input_path output_path')
    args = args_template(input_path=datadir[NOT_A_PAIR], output_path=outpath.strpath)
    main(args)
    with Reader(outpath.strpath, external_bin=None) as reader:
        assert len([r for r in reader if r.is_proper_pair]) == 2
    # Now do it again, but converting the namedtuple to an argv list
    argv = namedtuple_to_argv(args, 'allow_dovetailing.py')
    mocker.patch('sys.argv', argv)
    main()
    with Reader(outpath.strpath, external_bin=None) as reader:
        assert len([r for r in reader if r.is_proper_pair]) == 2
