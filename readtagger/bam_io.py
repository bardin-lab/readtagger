import gzip
import os
import subprocess
import pysam
import six
if six.PY2:
    import shutilwhich  # noqa: F401
import shutil  # noqa: E402


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
            if shutil.which('sambamba'):
                return sambamba_args
            if shutil.which('samtools'):
                return samtools_args
        return None

    def close(self):
        """Close filehandles and suprocess safely."""
        self.af.close()
        if self.args:
            self.proc.stdin.close()
            self.proc.wait()

    def __enter__(self):
        """Provide context handler entry."""
        if self.args:
            self.proc = subprocess.Popen(self.args, stdin=subprocess.PIPE, env=os.environ.copy(), close_fds=True)
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
            if shutil.which('sambamba'):
                return sambamba_args
            if shutil.which('samtools'):
                return samtools_args
        return None

    @property
    def is_bam(self):
        """Return whether file is BAM (True) or SAM."""
        if not hasattr(self, '_is_bam'):
            try:
                g = gzip.GzipFile(self.path)
                if g.read(3) == 'BAM':
                    self._is_bam = True
                else:
                    self._is_bam = False
            except Exception:
                self._is_bam = False
        return self._is_bam

    def close(self):
        """Close filehandles and suprocess safely."""
        self.af.close()
        if self.is_bam and self.args:
            self.proc.stdout.close()
            self.proc.wait()

    def __enter__(self):
        """Provide context handler entry."""
        if self.is_bam and self.args:
            self.proc = subprocess.Popen(self.args, stdout=subprocess.PIPE, env=os.environ.copy(), close_fds=True)
            self.af = pysam.AlignmentFile(self.proc.stdout)
        else:
            self.af = pysam.AlignmentFile(self.path)
        return self.af

    def __exit__(self, type, value, traceback):
        """Provide context handler exit."""
        self.close()
