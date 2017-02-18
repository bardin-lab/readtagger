from readtagger.allow_dovetailing import main
from readtagger.bam_io import BamAlignmentReader as Reader
from collections import namedtuple

NOT_A_PAIR = 'not_a_proper_pair.bam'


def test_main(datadir, tmpdir):  # noqa: D103
    outpath = tmpdir.join('pair.bam')
    args_template = namedtuple('ArgumentParser', 'input_path output_path')
    args = args_template(input_path=datadir[NOT_A_PAIR], output_path=outpath.strpath)
    main(args)
    with Reader(outpath.strpath, external_bin=None) as reader:
        assert len([r for r in reader if r.is_proper_pair]) == 2
