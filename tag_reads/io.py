import shutilwhich  # noqa: F401
from shutil import which
import subprocess
import pysam


class BamAlignmentWriter(object):
    """Wrap pysam.AlignmentFile with sambamba for multithreaded compressed writing."""

    def __init__(self, path, template=None, header=None, threads=2):
        """Use this class with a contexthandler."""
        self.template = template
        self.header = header
        self.path = path
        if which('sambamba'):
            self.args = ['sambamba', 'view', '-S', '-f', 'bam', '-t', "%s" % threads, '/dev/stdin', '-o', path]
        elif which('samtools'):
            self.args = ['samtools', 'view', '-b', '>', path]
        else:
            self.args = None

    def close(self):
        """Close filehandles and suprocess safely."""
        self.af.close()
        if self.args:
            self.proc.stdin.close()

    def __enter__(self):
        """Provide context handler entry."""
        if self.args:
            self.proc = subprocess.Popen(self.args, stdin=subprocess.PIPE)
            self.af = pysam.AlignmentFile(self.proc.stdin, mode="wh", template=self.template, header=self.header)
        else:
            self.af = pysam.AlignmentFile(self.path, mode="wb", template=self.template, header=self.header)
        return self.af

    def __exit__(self, type, value, traceback):
        """Provide context handler exit."""
        self.close()


class BamAlignmentReader(object):
    """Wraps pysam.AlignmentFile with sambamba for reading if input file is a bam file."""

    def __init__(self, path):
        """Use this class with a contexthandler."""
        self.bam = pysam.AlignmentFile(path).is_bam
        self.path = path
        if which('sambamba'):
            self.args = ['sambamba', 'view', '-h', path]
        elif which('samtools'):
            self.args = ['samtools', 'view', '-h', path]
        else:
            self.args = None

    def close(self):
        """Close filehandles and suprocess safely."""
        self.af.close()
        if self.bam and self.args:
            self.proc.stdout.close()

    def __enter__(self):
        """Provide context handler entry."""
        if self.bam and self.args:
            self.proc = subprocess.Popen(self.args, stdout=subprocess.PIPE)
            self.af = pysam.AlignmentFile(self.proc.stdout)
        else:
            self.af = pysam.AlignmentFile(self.path)
        return self.af

    def __exit__(self, type, value, traceback):
        """Provide context handler exit."""
        self.close()
