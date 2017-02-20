from collections import namedtuple
from readtagger.bam_io import BamAlignmentReader as Reader
from readtagger.mateoperations import main

TEST_BAM = 'dm6.bam'
TEMPLATE_ARGS = namedtuple('args', ['file_to_annotate', 'annotate_source', 'output_path', 'mate_sequence_tag'])


def test_mateoperations(datadir, tmpdir):  # noqa: D103
    out = tmpdir.join('out.bam')
    args = TEMPLATE_ARGS(file_to_annotate=datadir[TEST_BAM], annotate_source=datadir[TEST_BAM], output_path=out.strpath, mate_sequence_tag='MS')
    main(args)
    with Reader(out.strpath, external_bin=None) as reader:
        reads = [r for r in reader]
        assert len([True for r in reads if r.has_tag('MS')]) == 2
        assert reads[0].query_sequence == reads[1].get_tag('MS')
