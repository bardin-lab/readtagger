import subprocess
import time

from tag_reads.io import BamAlignmentReader
from tag_reads.io import BamAlignmentWriter

TEST_SAM = 'testsam_a.sam'
TEST_SAM_B = 'testsam_b.sam'


def test_bamreader_sambamba(datadir):  # noqa: D103
    _test_bamreader(datadir, 'sambamba')


def test_bamreader_samtools(datadir):  # noqa: D103
    _test_bamreader(datadir, 'samtools')


def test_bamreader_internal(datadir):  # noqa: D103
    _test_bamreader(datadir, None)


def test_bamwriter_sambamba(datadir, tmpdir):  # noqa: D103
    _test_bamwriter(datadir, tmpdir, 'sambamba')


def test_bamwriter_samtools(datadir, tmpdir):  # noqa: D103
    _test_bamwriter(datadir, tmpdir, 'samtools')


def test_bamwriter_internal(datadir, tmpdir):  # noqa: D103
    _test_bamwriter(datadir, tmpdir, None)


def _test_bamreader(datadir, external_bin):  # noqa: D103
    with BamAlignmentReader(datadir[TEST_SAM], external_bin) as reader:
        assert len([r for r in reader]) == 2


def _test_bamwriter(datadir, tmpdir, external_bin):  # noqa: D103
    outfile = tmpdir.join('out.bam')
    with BamAlignmentReader(datadir[TEST_SAM], external_bin) as reader, BamAlignmentWriter(outfile.strpath, template=reader, external_bin=external_bin) as out:
        for r in reader:
            out.write(r)
    time.sleep(1)  # Sleep a little to allow closing of handle
    out = subprocess.call(['samtools', 'view', outfile.strpath])
    assert out == 0