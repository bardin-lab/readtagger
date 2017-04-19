import gzip
import os
import subprocess
import tempfile
import pysam
import six

if six.PY2:
    import shutilwhich  # noqa: F401
import shutil  # noqa: E402


def is_file_coordinate_sorted(path):
    """Determine if first 10000 reads are coordinate sorted."""
    with pysam.AlignmentFile(path) as f:
        i = 0
        current_start = 0
        current_tid = 0
        for r in f:
            if i > 10000:
                return True
            if not r.is_unmapped:
                i += 1
                if r.reference_start >= current_start and r.tid >= current_tid:
                    current_start = r.reference_start
                    current_tid = r.tid
                    continue
                else:
                    return False
    return True


class BamAlignmentWriter(object):
    """Wrap pysam.AlignmentFile with sambamba for multithreaded compressed writing."""

    def __init__(self, path, template=None, header=None, threads=2, external_bin='choose_best', sort_order='coordinate'):
        """
        Write Bam files.

        Use this class with a context handler.

        :param path: Wite bam file to location `path`.
        :param template: Specify either template or header to write a bam file.
        :param header: Specify either template or header to write a bam file.
        :param external_bin: Specify `samtools` to use samtools or `sambamba` or `None` to use pysam for writing to `path`.
        :param sort_order: Can be `coordinate` or `queryname` and will cause the output file to sorted by this strategy.
        """
        self.template = template
        self.header = header
        self.path = path
        self.external_bin = external_bin
        self.threads = threads
        self.sort_order = sort_order

    @property
    def args(self):
        """Figure out if sambamba or samtools are available and return correct arguments."""
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
        """
        Close filehandles and subprocess safely.

        If necessary will sort the the file.
        """
        self.af.close()
        if self.args:
            self.proc.stdin.close()
            self.proc.wait()
        try:
            if is_file_coordinate_sorted(self.path):
                sort_order = 'coordinate'
            else:
                sort_order = 'queryname'
            if sort_order != self.sort_order:
                _, newpath = tempfile.mkstemp()
                args = ['samtools', 'sort']
                if self.sort_order == 'queryname':
                    args.append('-n')
                args.extend(['-o', newpath, self.path])
                subprocess.call(args, env=os.environ.copy())
                shutil.move(newpath, self.path)
        except Exception:
            # If no reads had been written to self.path
            # and exception will be reaised by is_file_coordinate_sorted.
            pass

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

    def __init__(self, path, external_bin='choose_best', sort_order='coordinate'):
        """
        Read Bam files.

        Use this class with a contexthandler.

        :param path: Path to read bam file from.
        :param external_bin: Specify `samtools` to use samtools or `sambamba` or `None` to use pysam for writing to `path`.
        :param sort_order: Can be `coordinate` or `queryname` and will cause the output file to sorted by this strategy.
        """
        self.path = path
        self.external_bin = external_bin
        self.sort_order = sort_order

    @property
    def args(self):
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
                if g.read(3) == 'BAM' or b'BAM':
                    self._is_bam = True
                else:
                    self._is_bam = False
            except Exception:
                self._is_bam = False
        return self._is_bam

    def close(self):
        """Close filehandles and subprocess safely."""
        self.af.close()
        if self.is_bam and self.args:
            self.proc.stdout.close()
            self.proc.wait()

    def __enter__(self):
        """Provide context handler entry."""
        if is_file_coordinate_sorted(self.path):
            sort_order = 'coordinate'
        else:
            sort_order = 'queryname'
        if sort_order != self.sort_order:
            _, newpath = tempfile.mkstemp()
            args = ['samtools', 'sort']
            if self.sort_order == 'queryname':
                args.append('-n')
            args.extend(['-o', newpath, self.path])
            subprocess.call(args, env=os.environ.copy())
            self.path = newpath
        if self.is_bam and self.args:
            self.proc = subprocess.Popen(self.args, stdout=subprocess.PIPE, env=os.environ.copy(), close_fds=True)
            self.af = pysam.AlignmentFile(self.proc.stdout)
        else:
            self.af = pysam.AlignmentFile(self.path)
        return self.af

    def __exit__(self, type, value, traceback):
        """Provide context handler exit."""
        self.close()
