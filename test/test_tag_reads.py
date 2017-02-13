import subprocess

from tag_reads.io import BamAlignmentReader
from tag_reads.io import BamAlignmentWriter

from tag_reads.tag_reads import SamTagProcessor
from tag_reads.tag_reads import SamAnnotator

TEST_SAM = 'testsam_a.sam'
TEST_SAM_B = 'testsam_b.sam'


def test_bamreader(datadir):
    with BamAlignmentReader(datadir[TEST_SAM]) as reader:
        assert len([r for r in reader]) == 2


def test_bamwriter(datadir, tmpdir):
    outfile = tmpdir.join('out.bam')
    with BamAlignmentReader(datadir[TEST_SAM]) as reader, BamAlignmentWriter(outfile.strpath, template=reader) as out:
        for r in reader:
            out.write(r)
    out = subprocess.call(['sambamba', 'view', outfile.strpath])
    assert out == 0


def test_samtag_processor(datadir):
    p_mate_false = get_samtag_processor(datadir, tag_mate=False)
    assert len(p_mate_false.result) == 1
    assert len(p_mate_false.result.values()) == 1
    assert True not in next(iter(p_mate_false.result.values()))
    p_mate_true = get_samtag_processor(datadir, tag_mate=True)
    assert next(iter(p_mate_true.result.values())) == next(iter(p_mate_false.result.values()))
    # Need some more data to make more meaningful tests


def test_samtag_annotator(datadir, tmpdir):
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


def get_samtag_processor(datadir, tag_mate):
    source_path = datadir[TEST_SAM]
    tag_prefix_self = 'A'
    tag_prefix_mate = 'B'
    tag_mate=False
    p = SamTagProcessor(source_path, tag_prefix_self, tag_prefix_mate, tag_mate)
    assert hasattr(p, 'result')
    assert isinstance(p.result, dict)
    return p