from readtagger.extract_variants import extract_variant_portion
import pysam

from .helpers import reference_fasta  # noqa: F401


def test_variant_portion(datadir_copy, tmpdir, reference_fasta):  # noqa: D103, F811
    input_path = str(datadir_copy['long_insertion.bam'])
    output_path = str(tmpdir.join('test_out.bam'))
    extract_variant_portion(input_path, output_path, reference_fasta)
    with pysam.AlignmentFile(output_path) as result:
        reads = [r for r in result]
        assert len(reads) == 1
        assert reads[0].cigarstring == '1M2317I1M'
        assert reads[0].query_sequence.startswith('N')
        assert reads[0].query_sequence.endswith('N')
