from readtagger.readtagger import SamTagProcessor
from readtagger.readtagger import SamAnnotator
from readtagger.readtagger import main
from .helpers import namedtuple_to_argv

from collections import namedtuple


TEST_SAM = 'testsam_a.sam'
TEST_SAM_B = 'testsam_b.sam'

TEST_BAM_A = 'dm6.bam'
TEST_BAM_B = 'pasteurianus.bam'
ARGS_TEMPLATE = namedtuple('args', ['annotate_with',
                                    'tag_file',
                                    'allow_dovetailing',
                                    'keep_suboptimal_alternate_tags',
                                    'discard_if_proper_pair',
                                    'write_discarded',
                                    'write_verified',
                                    'output_file'])


def test_samtag_processor(datadir):  # noqa: D103
    p_mate_false = get_samtag_processor(datadir, tag_mate=False)
    assert len(p_mate_false.result) == 1
    assert len(p_mate_false.result.values()) == 1
    assert True not in next(iter(p_mate_false.result.values()))
    p_mate_true = get_samtag_processor(datadir, tag_mate=True)
    assert next(iter(p_mate_true.result.values())) == next(iter(p_mate_false.result.values()))
    # Need some more data to make more meaningful tests


def test_samtag_annotator(datadir, tmpdir):  # noqa: D103
    p = get_samtag_processor(datadir, tag_mate=True)
    output_path = tmpdir.join('testout.bam')
    a = SamAnnotator(datadir[TEST_SAM_B],
                     [p],
                     output_path=output_path.strpath,
                     allow_dovetailing=False,
                     discard_bad_alt_tag=True,
                     discarded_writer=None,
                     verified_writer=None)
    assert isinstance(a, SamAnnotator)
    output_path.check()


def test_main_keep_suboptimal(datadir, tmpdir):  # noqa: D103
    # Annotate dm6 with pasteurianus reads, keep suboptimal tags
    discarded, verified, output = get_output_files(tmpdir)
    annotate_with = str(datadir[TEST_BAM_B])
    tag_file = str(datadir[TEST_BAM_A])
    args = ARGS_TEMPLATE(annotate_with=[annotate_with], tag_file=tag_file, allow_dovetailing=True, keep_suboptimal_alternate_tags=True,
                         discard_if_proper_pair=False, output_file=output.strpath, write_discarded=discarded.strpath, write_verified=verified.strpath)
    main(args)


def test_main_discard_suboptimal(datadir, tmpdir):  # noqa: D103
    # Annotate dm6 with pasteurianus reads, keep suboptimal tags
    discarded, verified, output = get_output_files(tmpdir)
    annotate_with = str(datadir[TEST_BAM_B])
    tag_file = str(datadir[TEST_BAM_A])
    args = ARGS_TEMPLATE(annotate_with=[annotate_with], tag_file=tag_file, allow_dovetailing=True, keep_suboptimal_alternate_tags=False,
                         discard_if_proper_pair=False, output_file=output.strpath, write_discarded=discarded.strpath, write_verified=verified.strpath)
    main(args)


def test_main_discard_suboptimal_discard_if_proper(datadir, tmpdir):  # noqa: D103
    # Annotate dm6 with pasteurianus reads, keep suboptimal tags, discard proper pairs
    discarded, verified, output = get_output_files(tmpdir)
    annotate_with = str(datadir[TEST_BAM_B])
    tag_file = str(datadir[TEST_BAM_A])
    args = ARGS_TEMPLATE(annotate_with=[annotate_with], tag_file=tag_file, allow_dovetailing=True, keep_suboptimal_alternate_tags=False,
                         discard_if_proper_pair=True, output_file=output.strpath, write_discarded=discarded.strpath, write_verified=verified.strpath)
    main(args)


def test_main_with_argparse(datadir, tmpdir, mocker):  # noqa: D103
    # Annotate dm6 with pasteurianus reads, keep suboptimal tags, discard proper pairs
    # and test that argparse argument parsing works as expected.
    discarded, verified, output = get_output_files(tmpdir)
    annotate_with = str(datadir[TEST_BAM_B])
    tag_file = str(datadir[TEST_BAM_A])
    args = ARGS_TEMPLATE(annotate_with=[annotate_with], tag_file=tag_file, allow_dovetailing=True, keep_suboptimal_alternate_tags=False,
                         discard_if_proper_pair=True, output_file=output.strpath, write_discarded=discarded.strpath, write_verified=verified.strpath)
    argv = namedtuple_to_argv(args, 'readtagger.py')
    mocker.patch('sys.argv', argv)
    main()


def get_samtag_processor(datadir, tag_mate):  # noqa: D103
    source_path = datadir[TEST_SAM]
    tag_prefix_self = 'A'
    tag_prefix_mate = 'B'
    p = SamTagProcessor(source_path, tag_prefix_self, tag_prefix_mate, tag_mate)
    assert hasattr(p, 'result')
    assert isinstance(p.result, dict)
    return p


def get_output_files(tmpdir):  # noqa: D103
    return tmpdir.join('discarded.bam'), tmpdir.join('verified.bam'), tmpdir.join('output.bam')
