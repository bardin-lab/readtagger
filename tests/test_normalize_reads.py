from readtagger.normalization import (
    get_read_sizes,
    split_fastq_files,
)

FASTQ_A = 'long_reads_a.fastq'
FASTQ_B = 'long_reads_b.fastq'


def test_split_fastq_files(datadir_copy, tmpdir):  # noqa: D103
    input_files = [str(datadir_copy[FASTQ_A]), str(datadir_copy[FASTQ_B])]
    output_paths = [tmpdir.join(FASTQ_A).strpath, tmpdir.join(FASTQ_B).strpath]
    split_fastq_files(input_files, output_paths)
    assert get_read_sizes(output_paths[0]) == get_read_sizes(output_paths[1])
