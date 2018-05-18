import pytest
from readtagger.create_multisample_vcf import (
    VCFMerger,
    NotSingleSampleVcfException
)

TEST_VCF_FILES = ['h3.vcf.gz', 'h6.vcf.gz', 'h9.vcf.gz']
TEST_VCF_INDICES = ['h3.vcf.gz.csi', 'h6.vcf.gz.csi', 'h9.vcf.gz.csi']

TEST_BCF_FILES = ['H6.bcf', 'H9.bcf']
TEST_BCF_INDICES = ['H6.bcf.csi', 'H9.bcf.csi']

MULTISAMPLE_VCF = ['multisample.vcf.gz']
MULTISAMPLE_VCF_INDEX = ['multisample.vcf.gz.csi']


def test_vcf_merger_vcf(datadir_copy, tmpdir):  # noqa: D103
    variant_files = [str(datadir_copy[f]) for f in TEST_VCF_FILES]
    [str(datadir_copy[f]) for f in TEST_VCF_INDICES]
    output_path = tmpdir.join('merged.vcf').strpath
    VCFMerger(variant_file_paths=variant_files, output_path=output_path, window_size=2)


def test_vcf_merger_bcf(datadir_copy, tmpdir):  # noqa: D103
    variant_files = [str(datadir_copy[f]) for f in TEST_VCF_FILES]
    [str(datadir_copy[f]) for f in TEST_VCF_INDICES]
    output_path = tmpdir.join('merged.vcf').strpath
    VCFMerger(variant_file_paths=variant_files, output_path=output_path, window_size=300, search_window=-500)


def test_vcf_merger_bcf2(datadir_copy, tmpdir):  # noqa: D103
    variant_files = [str(datadir_copy[f]) for f in TEST_BCF_FILES]
    [str(datadir_copy[f]) for f in TEST_BCF_INDICES]
    output_path = tmpdir.join('merged.vcf').strpath
    VCFMerger(variant_file_paths=variant_files, output_path=output_path, window_size=300, search_window=300)


def test_vcf_merger_multisample_exception(datadir_copy, tmpdir):  # noqa: D103
    variant_files = [str(datadir_copy[f]) for f in MULTISAMPLE_VCF]
    [str(datadir_copy[f]) for f in MULTISAMPLE_VCF_INDEX]
    output_path = tmpdir.join('merged.vcf').strpath
    with pytest.raises(NotSingleSampleVcfException):
        VCFMerger(variant_file_paths=variant_files, output_path=output_path, window_size=300, search_window=300)
