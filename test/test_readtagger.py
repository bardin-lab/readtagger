from readtagger.readtagger import SamTagProcessor
from readtagger.readtagger import SamAnnotator
from readtagger.readtagger import main
from readtagger.bam_io import (
    BamAlignmentReader as Reader,
    BamAlignmentWriter as Writer,
)
from .helpers import namedtuple_to_argv

import pysam
from collections import namedtuple


TEST_SAM = 'testsam_a.sam'
TEST_SAM_B = 'testsam_b.sam'

TEST_BAM_A = 'dm6.bam'
TEST_BAM_B = 'pasteurianus.bam'

TEST_SAM_ROVER_DM6 = 'rover_single_mate_dm6.sam'
TEST_SAM_ROVER_FBTI = 'rover_single_mate_fbti.sam'

ARGS_TEMPLATE = namedtuple('args', ['annotate_with',
                                    'tag_file',
                                    'allow_dovetailing',
                                    'keep_suboptimal_alternate_tags',
                                    'discard_if_proper_pair',
                                    'write_discarded',
                                    'write_verified',
                                    'output_file',
                                    'cores'])


def test_samtag_processor(datadir):  # noqa: D103
    # test that tag_mate=False really skips tagging mates
    p_mate_false = get_samtag_processor(datadir, tag_mate=False)
    for qname, tag_d in p_mate_false.process():
        assert qname == 'DISCARD_R1_KEEP_R2:1'
        assert len(tag_d) == 1
    p_mate_true = get_samtag_processor(datadir, tag_mate=True)
    # Make sure that the mate gets tagged as well
    for qname, tag_d in p_mate_true.process():
        assert len(tag_d) == 2
    # Need some more data to make more meaningful tests


def test_samtag_annotator(datadir, tmpdir):  # noqa: D103
    p = get_samtag_processor(datadir, tag_mate=True)
    output_path = tmpdir.join('testout.bam')
    with Reader(datadir[TEST_SAM_B], sort_order='queryname') as annotate_bam:
        header = annotate_bam.header
        annotate_reads = [r for r in annotate_bam]
    with Writer(output_path.strpath, header=header) as output_writer:
        a = SamAnnotator(annotate_bam=annotate_reads,
                         output_writer=output_writer,
                         allow_dovetailing=False,
                         discard_bad_alt_tag=True,)
        assert isinstance(a, SamAnnotator)
        for qname, tag_d in p.process():
            a.process(qname, tag_d)
    output_path.check()


def test_main_keep_suboptimal(datadir, tmpdir):  # noqa: D103
    # Annotate dm6 with pasteurianus reads, keep suboptimal tags
    discarded, verified, output = get_output_files(tmpdir)
    annotate_with = str(datadir[TEST_BAM_B])
    tag_file = str(datadir[TEST_BAM_A])
    args = ARGS_TEMPLATE(annotate_with=[annotate_with], tag_file=tag_file, allow_dovetailing=True, keep_suboptimal_alternate_tags=True,
                         discard_if_proper_pair=False, output_file=output.strpath, write_discarded=discarded.strpath, write_verified=verified.strpath,
                         cores=1)
    main(args)


def test_main_discard_suboptimal(datadir, tmpdir):  # noqa: D103
    # Annotate dm6 with pasteurianus reads, keep suboptimal tags
    discarded, verified, output = get_output_files(tmpdir)
    annotate_with = str(datadir[TEST_BAM_B])
    tag_file = str(datadir[TEST_BAM_A])
    args = ARGS_TEMPLATE(annotate_with=[annotate_with], tag_file=tag_file, allow_dovetailing=True, keep_suboptimal_alternate_tags=False,
                         discard_if_proper_pair=False, output_file=output.strpath, write_discarded=discarded.strpath, write_verified=verified.strpath,
                         cores=1)
    main(args)


def test_main_discard_suboptimal_discard_if_proper(datadir, tmpdir):  # noqa: D103
    # Annotate dm6 with pasteurianus reads, keep suboptimal tags, discard proper pairs
    discarded, verified, output = get_output_files(tmpdir)
    annotate_with = str(datadir[TEST_BAM_B])
    tag_file = str(datadir[TEST_BAM_A])
    args = ARGS_TEMPLATE(annotate_with=[annotate_with], tag_file=tag_file, allow_dovetailing=True, keep_suboptimal_alternate_tags=False,
                         discard_if_proper_pair=True, output_file=output.strpath, write_discarded=discarded.strpath, write_verified=verified.strpath,
                         cores=1)
    main(args)


def test_main_with_argparse(datadir, tmpdir, mocker):  # noqa: D103
    # Annotate dm6 with pasteurianus reads, keep suboptimal tags, discard proper pairs
    # and test that argparse argument parsing works as expected.
    discarded, verified, output = get_output_files(tmpdir)
    annotate_with = str(datadir[TEST_BAM_B])
    tag_file = str(datadir[TEST_BAM_A])
    args = ARGS_TEMPLATE(annotate_with=[annotate_with], tag_file=tag_file, allow_dovetailing=True, keep_suboptimal_alternate_tags=False,
                         discard_if_proper_pair=True, output_file=output.strpath, write_discarded=discarded.strpath, write_verified=verified.strpath,
                         cores='1')
    argv = namedtuple_to_argv(args, 'readtagger.py')
    mocker.patch('sys.argv', argv)
    main()


def test_main_rover(datadir, tmpdir, mocker):  # noqa: D103
    # Annotate dm6 with pasteurianus reads, keep suboptimal tags, discard proper pairs
    # and test that argparse argument parsing works as expected.
    discarded, verified, output = get_output_files(tmpdir)
    annotate_with = str(datadir[TEST_SAM_ROVER_FBTI])
    tag_file = str(datadir[TEST_SAM_ROVER_DM6])
    args = ARGS_TEMPLATE(annotate_with=[annotate_with], tag_file=tag_file, allow_dovetailing=True, keep_suboptimal_alternate_tags=False,
                         discard_if_proper_pair=True, output_file=output.strpath, write_discarded=discarded.strpath, write_verified=verified.strpath,
                         cores='1')
    argv = namedtuple_to_argv(args, 'readtagger.py')
    mocker.patch('sys.argv', argv)
    main()
    assert len([r for r in pysam.AlignmentFile(verified.strpath)]) == 1


def get_samtag_processor(datadir, tag_mate):  # noqa: D103
    source_path = datadir[TEST_SAM]
    header = pysam.AlignmentFile(source_path).header
    return SamTagProcessor(Reader(source_path, sort_order='queryname').__enter__(), header=header, tag_mate=tag_mate)


def get_output_files(tmpdir):  # noqa: D103
    return tmpdir.join('discarded.bam'), tmpdir.join('verified.bam'), tmpdir.join('output.bam')
