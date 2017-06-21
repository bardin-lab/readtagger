import pysam

from readtagger.pysamtools_view import view

INPUT = 'tagged_dm6.bam'


def test_pysamtoolsview(datadir, tmpdir):  # noqa: D103
    input_bam = datadir[INPUT]
    output_bam = tmpdir.join('out.bam').strpath
    region = '2L:2878990-2880990'
    view(input_bam=input_bam, output_bam=output_bam, region=region)
    assert len(pysam.AlignmentFile(output_bam).header['SQ']) == 1
