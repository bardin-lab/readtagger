from readtagger.bam_io import BamAlignmentReader
from readtagger.bam_io import BamAlignmentWriter

TEST_BAM = 'dm6.bam'
TEST_SAM = 'testsam_a.sam'
TEST_SAM_B = 'testsam_b.sam'


def test_bamreader_sambamba(datadir):  # noqa: D103
    _test_bamreader(datadir, 'sambamba')


def test_bamreader_samtools(datadir):  # noqa: D103
    _test_bamreader(datadir, 'samtools')


def test_bamreader_internal(datadir):  # noqa: D103
    _test_bamreader(datadir, None)


def test_bamreader_choose_best_samtools(datadir, mocker):  # noqa: D103
    mocker.patch('shutil.which', side_effect=return_samtools)  # will test samtools
    _test_bamreader(datadir, 'choose_best')


def test_bamreader_bam(datadir):  # noqa: D103
    with BamAlignmentReader(datadir[TEST_BAM], 'sambamba') as reader:
        assert len([r for r in reader]) == 2


def test_bamwriter_sambamba(datadir, tmpdir):  # noqa: D103
    for i in range(10):  # Bit of stress testing here.
        _test_bamwriter(datadir, tmpdir, 'sambamba', i)


def test_bamwriter_samtools(datadir, tmpdir):  # noqa: D103
    for i in range(10):  # Bit of stress testing here.
        _test_bamwriter(datadir, tmpdir, 'samtools', i)


def test_bamwriter_internal(datadir, tmpdir):  # noqa: D103
    _test_bamwriter(datadir, tmpdir, None)


def test_bamwriter_choose_best_samtools(datadir, tmpdir, mocker):  # noqa: D103
    mocker.patch('shutil.which', side_effect=return_samtools)  # will test samtools
    _test_bamwriter(datadir, tmpdir, 'choose_best')


def _test_bamreader(datadir, external_bin):  # noqa: D103
    with BamAlignmentReader(datadir[TEST_SAM], external_bin) as reader:
        assert len([r for r in reader]) == 2


def _test_bamwriter(datadir, tmpdir, external_bin='choose_best', i=0):  # noqa: D103
    outfile = tmpdir.join("%s_%s" % (i, 'out.bam'))
    with BamAlignmentReader(datadir[TEST_SAM], external_bin) as reader, BamAlignmentWriter(outfile.strpath, template=reader, external_bin=external_bin) as out:
        for r in reader:
            out.write(r)
    with BamAlignmentReader(outfile.strpath, None) as reader:
        assert len([r for r in reader]) == 2


def return_samtools(arg):  # noqa: D103
    if arg == 'samtools':
        return True
    else:
        return False
