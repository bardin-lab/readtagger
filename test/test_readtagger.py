from readtagger.readtagger import SamTagProcessor
from readtagger.readtagger import SamAnnotator
from readtagger.readtagger import main

from collections import namedtuple


TEST_SAM = 'testsam_a.sam'
TEST_SAM_B = 'testsam_b.sam'

TEST_BAM_A = 'dm6.bam'
TEST_BAM_B = 'pasteurianus.bam'


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


def test_main(datadir, tmpdir):  # noqa: D103
    discarded = tmpdir.join('discarded.bam')
    verified = tmpdir.join('verified.bam')
    output = tmpdir.join('output.bam')
    annotate_with = str(datadir[TEST_BAM_A])
    tag_file = str(datadir[TEST_BAM_B])
    args_template = namedtuple('args', 'annotate_with tag_file allow_dovetailing keep_suboptimal_alternate_tags write_discarded write_verified, output_file')
    args = args_template(annotate_with=[annotate_with], tag_file=tag_file, allow_dovetailing=True, keep_suboptimal_alternate_tags=True,
                         output_file=output.strpath, write_discarded=discarded.strpath, write_verified=verified.strpath)
    main(args)


def get_samtag_processor(datadir, tag_mate):  # noqa: D103
    source_path = datadir[TEST_SAM]
    tag_prefix_self = 'A'
    tag_prefix_mate = 'B'
    p = SamTagProcessor(source_path, tag_prefix_self, tag_prefix_mate, tag_mate)
    assert hasattr(p, 'result')
    assert isinstance(p.result, dict)
    return p
