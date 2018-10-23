import pysam
import readtagger.bam_io


TEST_BAM = 'dm6.bam'
TEST_SAM = 'testsam_a.sam'
TEST_SAM_B = 'testsam_b.sam'
EXTENDED = 'extended_annotated_updated_all_reads.bam'


def test_bamreader_samtools(datadir_copy):  # noqa: D103
    _test_bamreader(datadir_copy, 'samtools')


def test_bamreader_internal(datadir_copy):  # noqa: D103
    _test_bamreader(datadir_copy, None)


def test_bamreader_queryname_sort(datadir_copy):  # noqa: D103
    _test_bamreader(datadir_copy, sort_order='queryname')


def test_bamreader_bam(datadir_copy):  # noqa: D103
    with readtagger.bam_io.BamAlignmentReader(str(datadir_copy[TEST_BAM]), 'samtools') as reader:
        assert len([r for r in reader]) == 2


def test_bamreader_region(datadir_copy):  # noqa: D103
    with readtagger.bam_io.BamAlignmentReader(str(datadir_copy[TEST_BAM]), 'samtools', region='2L:2877890-2878090') as reader:
        assert len([r for r in reader]) == 2


def test_bamreader_sam_threads(datadir_copy):  # noqa: D103
    with readtagger.bam_io.BamAlignmentReader(str(datadir_copy[TEST_SAM]), external_bin='samtools', threads=1) as reader:
        assert len([r for r in reader]) == 2


def test_bamwriter_samtools(datadir_copy, tmpdir):  # noqa: D103
    for i in range(10):  # Bit of stress testing here.
        _test_bamwriter(datadir_copy, tmpdir, 'samtools', i)


def test_bamwriter_internal(datadir_copy, tmpdir):  # noqa: D103
    _test_bamwriter(datadir_copy, tmpdir, None)


def test_bamwriter_queryname_sort(datadir_copy, tmpdir):  # noqa: D103
    _test_bamwriter(datadir_copy, tmpdir, sort_order='queryname')


def _test_bamreader(datadir_copy, external_bin='samtools', sort_order='coordinate'):  # noqa: D103
    with readtagger.bam_io.BamAlignmentReader(str(datadir_copy[TEST_SAM]), external_bin, sort_order=sort_order) as reader:
        assert len([r for r in reader]) == 2


def _test_bamwriter(datadir_copy, tmpdir, i=0, external_bin='samtools', sort_order='coordinate', threads=1):  # noqa: D103
    outfile = tmpdir.join("%s_%s" % (i, 'out.bam'))
    with readtagger.bam_io.BamAlignmentReader(str(datadir_copy[TEST_SAM]), external_bin, sort_order=sort_order) as reader, \
            readtagger.bam_io.BamAlignmentWriter(outfile.strpath,
                                                 template=reader,
                                                 sort_order=sort_order,
                                                 threads=threads) as out:
        for r in reader:
            out.write(r)
    with readtagger.bam_io.BamAlignmentReader(outfile.strpath, None) as reader:
        assert len([r for r in reader]) == 2


def test_bamwriter_empty(datadir_copy, tmpdir):  # noqa: D103
    outfile = tmpdir.join('out.bam')
    with readtagger.bam_io.BamAlignmentReader(str(datadir_copy[TEST_SAM]), external_bin='samtools', sort_order='coordinate') as reader, \
            readtagger.bam_io.BamAlignmentWriter(outfile.strpath, template=reader, sort_order='queryname') as _:  # noqa: F841
        pass


def test_sort_cram(datadir_copy, tmpdir):  # noqa: D103
    outfile = tmpdir.join('out.cram').strpath
    in_path = str(datadir_copy[EXTENDED])
    readtagger.bam_io.sort_bam(inpath=in_path, output=outfile, sort_order='queryname', cram=True)
    assert len([r for r in pysam.AlignmentFile(in_path)]) == len([r for r in pysam.AlignmentFile(outfile)])


def test_merge_bam(datadir_copy, tmpdir):  # noqa: D103
    in_path = str(datadir_copy[EXTENDED])
    outfile = tmpdir.join('out.cram').strpath
    bam_collection = [in_path, in_path]
    readtagger.bam_io.merge_bam(bam_collection, output_path=outfile)


def test_bamwriter_switch_output_sorting(datadir_copy, tmpdir):  # noqa: D103
    outfile = tmpdir.join('out.bam')
    with readtagger.bam_io.BamAlignmentReader(str(datadir_copy[EXTENDED]), external_bin='samtools', sort_order='coordinate', threads=1) as reader, \
            readtagger.bam_io.BamAlignmentWriter(outfile.strpath, template=reader, sort_order='queryname', threads=1) as writer:
        for r in reader:
            writer.write(r)
    with readtagger.bam_io.BamAlignmentReader(str(datadir_copy[EXTENDED]), external_bin='samtools', sort_order='queryname', threads=5) as reader, \
            readtagger.bam_io.BamAlignmentWriter(outfile.strpath, template=reader, sort_order='coordinate', threads=5) as writer:
        for r in reader:
            writer.write(r)


def test_is_file_coordinate_sorted(datadir_copy, tmpdir):  # noqa: D103
    qname_sorted = tmpdir.join('qname_sorted').strpath
    qname_sorted = readtagger.bam_io.sort_bam(inpath=str(datadir_copy[EXTENDED]), output=qname_sorted, sort_order='queryname')
    assert not readtagger.bam_io.is_file_coordinate_sorted(qname_sorted, reads_to_check=10)


def test_get_queryname_positions(datadir_copy, tmpdir):  # noqa: D103
    qname_sorted = tmpdir.join('qname_sorted').strpath
    qname_sorted = readtagger.bam_io.sort_bam(inpath=str(datadir_copy[EXTENDED]), output=qname_sorted, sort_order='queryname')
    assert readtagger.bam_io.is_file_coordinate_sorted(str(datadir_copy[EXTENDED]), reads_to_check=10)
    assert not readtagger.bam_io.is_file_coordinate_sorted(qname_sorted, reads_to_check=10)
    qname_pos = readtagger.bam_io.get_queryname_positions(qname_sorted, chunk_size=50)
    assert len(qname_pos) == 8
    start_positions = [t[0] for t in qname_pos]
    last_qnames = [t[1] for t in qname_pos]
    start_positions_for_last_qnames = readtagger.bam_io.start_positions_for_last_qnames(qname_sorted, last_qnames[::])
    assert start_positions_for_last_qnames == start_positions
    readtagger.bam_io.get_reads(qname_sorted, start=start_positions[0], last_qname=last_qnames[0])
    readtagger.bam_io.get_reads(qname_sorted, start=start_positions[1], last_qname=last_qnames[1])
    start_positions_for_last_qnames = readtagger.bam_io.start_positions_for_last_qnames(qname_sorted, last_qnames[:1])
    assert len(start_positions_for_last_qnames) == 2


def test_start_positions_for_last_qnames(datadir_copy):  # noqa: D103
    bam = str(datadir_copy[EXTENDED])
    r = readtagger.bam_io.start_positions_for_last_qnames(bam, ['i_dont_exist'])
    assert len(r) == 1  # that's the start position


def test_find_end(datadir_copy):  # noqa: D103
    bam = str(datadir_copy[EXTENDED])
    with readtagger.bam_io.BamAlignmentReader(bam, external_bin=False, index=True) as f:
        end = readtagger.bam_io.find_end(f=f, chrom='3R', self_tag='AD', other_tag='BD', end=13373438, padding=100)
        assert end == 13373366


def test_split_locations_between_clusters(datadir_copy, tmpdir):  # noqa: D103
    bam = str(datadir_copy[EXTENDED])
    region = '3R:13372696-13373943'
    r = readtagger.bam_io.split_locations_between_clusters(bam, self_tag='AD', other_tag='BD', distance=100, region=region)
    assert len(r) == 9


def test_get_mean_read_length(datadir_copy):  # noqa: D103
    mean_rl = readtagger.bam_io.get_mean_read_length(str(datadir_copy[EXTENDED]), reads_to_check=10)
    assert mean_rl == 125


def return_samtools(arg):  # noqa: D103
    if arg == 'samtools':
        return True
    else:
        return False
