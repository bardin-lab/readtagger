import pysam
from readtagger.vcf_io import merge_vcf_files

TEST_VCF_FILE = 'test_vcf.vcf'


def test_merge_vcf_files(datadir_copy, tmpdir):  # noqa: D103
    test_vcf = str(datadir_copy[TEST_VCF_FILE])
    output_vcf = tmpdir.join('output.vcf').strpath
    merge_vcf_files([test_vcf, test_vcf], output_path=output_vcf)
    with pysam.VariantFile(test_vcf) as vf:
        original_count = len(list(vf))
    with pysam.VariantFile(output_vcf) as vf:
        merged_count = len(list(vf))
    assert merged_count == 2 * original_count
