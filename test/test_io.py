import readtagger.bam_io


TEST_BAM = 'dm6.bam'
TEST_SAM = 'testsam_a.sam'
TEST_SAM_B = 'testsam_b.sam'
EXTENDED = 'extended_annotated_updated_all_reads.bam'


def test_bamreader_samtools(datadir):  # noqa: D103
    _test_bamreader(datadir, 'samtools')


def test_bamreader_internal(datadir):  # noqa: D103
    _test_bamreader(datadir, None)


def test_bamreader_queryname_sort(datadir):  # noqa: D103
    _test_bamreader(datadir, sort_order='queryname')


def test_bamreader_bam(datadir):  # noqa: D103
    with readtagger.bam_io.BamAlignmentReader(datadir[TEST_BAM], 'samtools') as reader:
        assert len([r for r in reader]) == 2


def test_bamreader_region(datadir):  # noqa: D103
    with readtagger.bam_io.BamAlignmentReader(datadir[TEST_BAM], 'samtools', region='2L:2877890-2878090') as reader:
        assert len([r for r in reader]) == 2


def test_bamreader_sam_threads(datadir):  # noqa: D103
    with readtagger.bam_io.BamAlignmentReader(datadir[TEST_SAM], external_bin='samtools', threads=1) as reader:
        assert len([r for r in reader]) == 2


def test_bamwriter_samtools(datadir, tmpdir):  # noqa: D103
    for i in range(10):  # Bit of stress testing here.
        _test_bamwriter(datadir, tmpdir, 'samtools', i)


def test_bamwriter_internal(datadir, tmpdir):  # noqa: D103
    _test_bamwriter(datadir, tmpdir, None)


def test_bamwriter_queryname_sort(datadir, tmpdir):  # noqa: D103
    _test_bamwriter(datadir, tmpdir, sort_order='query_name')


def _test_bamreader(datadir, external_bin='samtools', sort_order='coordinate'):  # noqa: D103
    with readtagger.bam_io.BamAlignmentReader(datadir[TEST_SAM], external_bin, sort_order=sort_order) as reader:
        assert len([r for r in reader]) == 2


def _test_bamwriter(datadir, tmpdir, i=0, external_bin='samtools', sort_order='coordinate'):  # noqa: D103
    outfile = tmpdir.join("%s_%s" % (i, 'out.bam'))
    with readtagger.bam_io.BamAlignmentReader(datadir[TEST_SAM], external_bin, sort_order=sort_order) as reader, \
            readtagger.bam_io.BamAlignmentWriter(outfile.strpath, template=reader, external_bin=external_bin, sort_order=sort_order) as out:
        for r in reader:
            out.write(r)
    with readtagger.bam_io.BamAlignmentReader(outfile.strpath, None) as reader:
        assert len([r for r in reader]) == 2


def test_bamwriter_empty(datadir, tmpdir):  # noqa: D103
    outfile = tmpdir.join('out.bam')
    with readtagger.bam_io.BamAlignmentReader(datadir[TEST_SAM], external_bin='samtools', sort_order='coordinate') as reader, \
            readtagger.bam_io.BamAlignmentWriter(outfile.strpath, template=reader, external_bin='samtools', sort_order='queryname') as _:  # noqa: F841
        pass


def test_bamwriter_switch_output_sorting(datadir, tmpdir):  # noqa: D103
    outfile = tmpdir.join('out.bam')
    with readtagger.bam_io.BamAlignmentReader(datadir[EXTENDED], external_bin='samtools', sort_order='coordinate', threads=1) as reader, \
            readtagger.bam_io.BamAlignmentWriter(outfile.strpath, template=reader, external_bin='samtools', sort_order='queryname', threads=1) as writer:
        for r in reader:
            writer.write(r)
    with readtagger.bam_io.BamAlignmentReader(datadir[EXTENDED], external_bin='samtools', sort_order='queryname', threads=5) as reader, \
            readtagger.bam_io.BamAlignmentWriter(outfile.strpath, template=reader, external_bin='samtools', sort_order='coordinate', threads=5) as writer:
        for r in reader:
            writer.write(r)


def test_get_queryname_positions(datadir, tmpdir):  # noqa: D103
    qname_sorted = tmpdir.join('qname_sorted').strpath
    qname_sorted = readtagger.bam_io.sort_bam(inpath=datadir[EXTENDED], output=qname_sorted, sort_order='queryname')
    assert readtagger.bam_io.is_file_coordinate_sorted(datadir[EXTENDED], reads_to_check=10)
    assert not readtagger.bam_io.is_file_coordinate_sorted(qname_sorted, reads_to_check=10)
    qname_pos = readtagger.bam_io.get_queryname_positions(qname_sorted, chunk_size=50)
    assert len(qname_pos) == 8
    start_positions = [t[0] for t in qname_pos]
    last_qnames = [t[1] for t in qname_pos]
    start_positions_for_last_qnames = readtagger.bam_io.start_positions_for_last_qnames(qname_sorted, last_qnames[::])
    assert start_positions_for_last_qnames == start_positions
    readtagger.bam_io.get_reads(qname_sorted, start=start_positions[0], last_qname=last_qnames[0])
    readtagger.bam_io.get_reads(qname_sorted, start=start_positions[1], last_qname=last_qnames[1])


def return_samtools(arg):  # noqa: D103
    if arg == 'samtools':
        return True
    else:
        return False
