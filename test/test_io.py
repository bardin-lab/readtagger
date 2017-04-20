from readtagger.bam_io import BamAlignmentReader
from readtagger.bam_io import BamAlignmentWriter

TEST_BAM = 'dm6.bam'
TEST_SAM = 'testsam_a.sam'
TEST_SAM_B = 'testsam_b.sam'


def test_bamreader_samtools(datadir):  # noqa: D103
    _test_bamreader(datadir, 'samtools')


def test_bamreader_internal(datadir):  # noqa: D103
    _test_bamreader(datadir, None)


def test_bamreader_queryname_sort(datadir):  # noqa: D103
    _test_bamreader(datadir, sort_order='queryname')


def test_bamreader_bam(datadir):  # noqa: D103
    with BamAlignmentReader(datadir[TEST_BAM], 'samtools') as reader:
        assert len([r for r in reader]) == 2


def test_bamwriter_samtools(datadir, tmpdir):  # noqa: D103
    for i in range(10):  # Bit of stress testing here.
        _test_bamwriter(datadir, tmpdir, 'samtools', i)


def test_bamwriter_internal(datadir, tmpdir):  # noqa: D103
    _test_bamwriter(datadir, tmpdir, None)


def test_bamwriter_queryname_sort(datadir, tmpdir):  # noqa: D103
    _test_bamwriter(datadir, tmpdir, sort_order='query_name')


def _test_bamreader(datadir, external_bin='samtools', sort_order='coordinate'):  # noqa: D103
    with BamAlignmentReader(datadir[TEST_SAM], external_bin, sort_order=sort_order) as reader:
        assert len([r for r in reader]) == 2


def _test_bamwriter(datadir, tmpdir, i=0, external_bin='samtools', sort_order='coordinate'):  # noqa: D103
    outfile = tmpdir.join("%s_%s" % (i, 'out.bam'))
    with BamAlignmentReader(datadir[TEST_SAM], external_bin, sort_order=sort_order) as reader, \
            BamAlignmentWriter(outfile.strpath, template=reader, external_bin=external_bin, sort_order=sort_order) as out:
        for r in reader:
            out.write(r)
    with BamAlignmentReader(outfile.strpath, None) as reader:
        assert len([r for r in reader]) == 2


def return_samtools(arg):  # noqa: D103
    if arg == 'samtools':
        return True
    else:
        return False
