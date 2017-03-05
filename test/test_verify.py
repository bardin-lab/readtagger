from readtagger.verify import discard_supplementary
from readtagger.bam_io import BamAlignmentReader as Reader

SUP = 'supplementary.sam'


def test_verify(tmpdir, datadir):  # noqa: D103
    output_path = tmpdir.join('out.bam').strpath
    input_path = datadir[SUP]
    discard_supplementary(input_path=input_path, output_path=output_path)
    with Reader(output_path) as reader:
        assert len([r for r in reader]) == 0
