import shutilwhich  # noqa: F401
from shutil import which
import os
import subprocess
import pysam


class BamAlignmentWriter(object):
    """Wrap pysam.AlignmentFile with sambamba for multithreaded compressed writing."""

    def __init__(self, path, template=None, header=None, threads=2, external_bin='choose_best'):
        """Use this class with a contexthandler."""
        self.template = template
        self.header = header
        self.path = path
        self.external_bin = external_bin
        self.threads = threads
        self.args = self.get_subprocess_args()

    def get_subprocess_args(self):
        """Figure out in sambamba or samtools are available and return correct arguments."""
        sambamba_args = ['sambamba', 'view', '-f', 'bam', '-t', "%s" % self.threads, '/dev/stdin', '-o', self.path]
        samtools_args = ['samtools', 'view', '-b', '/dev/stdin', '-o', self.path]
        if self.external_bin == 'sambamba':
            return sambamba_args
        elif self.external_bin == 'samtools':
            return samtools_args
        elif self.external_bin == 'choose_best':
            if which('sambamba'):
                return sambamba_args
            if which('samtools'):
                return samtools_args
        return None

    def close(self):
        """Close filehandles and suprocess safely."""
        self.af.close()
        if self.args:
            self.proc.stdin.close()

    def __enter__(self):
        """Provide context handler entry."""
        if self.args:
            self.proc = subprocess.Popen(self.args, stdin=subprocess.PIPE, env=os.environ.copy())
            self.af = pysam.AlignmentFile(self.proc.stdin, mode="wbu", template=self.template, header=self.header)
        else:
            self.af = pysam.AlignmentFile(self.path, mode="wb", template=self.template, header=self.header)
        return self.af

    def __exit__(self, type, value, traceback):
        """Provide context handler exit."""
        self.close()


class BamAlignmentReader(object):
    """Wraps pysam.AlignmentFile with sambamba for reading if input file is a bam file."""

    def __init__(self, path, external_bin='choose_best'):
        """Use this class with a contexthandler."""
        self.bam = pysam.AlignmentFile(path).is_bam
        self.path = path
        self.external_bin = external_bin
        self.args = self.get_subprocess_args()

    def get_subprocess_args(self):
        """Figure out in sambamba or samtools are available and return correct arguments."""
        sambamba_args = ['sambamba', 'view', '-h', self.path]
        samtools_args = ['samtools', 'view', '-h', self.path]
        if self.external_bin == 'sambamba':
            return sambamba_args
        elif self.external_bin == 'samtools':
            return samtools_args
        elif self.external_bin == 'choose_best':
            if which('sambamba'):
                return sambamba_args
            if which('samtools'):
                return samtools_args
        return None

    def close(self):
        """Close filehandles and suprocess safely."""
        self.af.close()
        if self.bam and self.args:
            self.proc.stdout.close()

    def __enter__(self):
        """Provide context handler entry."""
        if self.bam and self.args:
            self.proc = subprocess.Popen(self.args, stdout=subprocess.PIPE, env=os.environ.copy())
            self.af = pysam.AlignmentFile(self.proc.stdout)
        else:
            self.af = pysam.AlignmentFile(self.path)
        return self.af

    def __exit__(self, type, value, traceback):
        """Provide context handler exit."""
        self.close()
