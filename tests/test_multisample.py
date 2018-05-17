from readtagger.create_multisample_vcf import VCFMerger

TEST_VCF_FILES = ['h3.vcf.gz', 'h6.vcf.gz', 'h9.vcf.gz']
TEST_VCF_INDICES = ['h3.vcf.gz.csi', 'h6.vcf.gz.csi', 'h9.vcf.gz.csi']


def test_vcf_merger(datadir_copy, tmpdir):  # noqa: D103
    variant_files = [str(datadir_copy[f]) for f in TEST_VCF_FILES]
    [str(datadir_copy[f]) for f in TEST_VCF_INDICES]
    output_path = tmpdir.join('merged.vcf').strpath
    VCFMerger(variant_file_paths=variant_files, output_path=output_path, window_size=2)
