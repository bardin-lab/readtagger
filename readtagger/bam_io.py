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


def get_queryname_positions(fn, chunk_size=10000):
    """
    Get positions in order to split alignment file into equal-sized portions.

    Return a list of tuples, where the first item is the current file position for the first read,
    and the second item is the last query_name of the chunk.
    """
    f = pysam.AlignmentFile(fn)
    start = f.tell()
    last_pos = start
    seek_positions = []
    r = next(f)
    qn = r.query_name
    count = 1
    for r in f:
        current_query_name = r.query_name
        if current_query_name != qn:
            # Next query_name
            count += 1
            if count % chunk_size == 0:
                # We've reach a chunk, we append the last_pos as a new start
                seek_positions.append((start, qn))
                start = last_pos
            qn = current_query_name
        last_pos = f.tell()
    seek_positions.append((start, qn))
    return seek_positions


def get_reads(fn, start, last_qname):
    """Get reads starting at `start` and ending with last_qname."""
    f = pysam.AlignmentFile(fn)
    reads = []
    f.seek(start)
    for r in f:
        if r.query_name == last_qname:
            try:
                while r.query_name == last_qname:
                    reads.append(r)
                    r = next(f)
            except StopIteration:
                pass
            return reads
        else:
            reads.append(r)
    return reads


def start_positions_for_last_qnames(fn, last_qnames):
    """Return start positions that return the first read after the current last qname."""
    f = pysam.AlignmentFile(fn)
    start = f.tell()
    current_last_qname = last_qnames.pop(0)
    seek_positions = [start]
    for r in f:
        if r.query_name == current_last_qname:
            last_file_pos = f.tell()
            try:
                while r.query_name == current_last_qname:
                    # We've got the last query_name,
                    # now we want the last file position in which current_last_qname occurs
                    last_file_pos = f.tell()
                    r = next(f)
                seek_positions.append(last_file_pos)
            except StopIteration:
                return seek_positions
            try:
                current_last_qname = last_qnames.pop(0)
            except IndexError:
                # We've reached the end of last_qnames
                return seek_positions
    return seek_positions


def merge_bam(bam_collection, template_bam, output_path, threads=1):
    """Merge a readname sorted collection of BAM files."""
    args = ['samtools', 'merge', '-n', '-f', '-@', "%s" % threads, '-h', template_bam, output_path]
    args.extend(bam_collection)
    subprocess.call(args, env=os.environ.copy())
    [os.remove(bam) for bam in bam_collection]
    return output_path


def sort_bam(inpath, output, sort_order, threads):
    """Sort bam file at inpath using sort_order and write output to output."""
    if inpath == output:
        # We sort 'in place'
        _, temp_out = tempfile.mkstemp()
    else:
        temp_out = output
    args = ['samtools', 'sort', '-@', "%s" % threads]
    if sort_order == 'queryname':
        args.append('-n')
    args.extend(['-o', temp_out, inpath])
    subprocess.call(args, env=os.environ.copy())
    if temp_out != output:
        shutil.move(temp_out, output)
    return output


class BamAlignmentWriter(object):
    """Wrap pysam.AlignmentFile with sambamba for multithreaded compressed writing."""

    def __init__(self, path, template=None, header=None, threads=4, external_bin='choose_best', sort_order='coordinate'):
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
        """Return samtools arguments for writing bam files."""
        if self.external_bin:
            if self.threads > 4:
                threads = 4  # More threads won't speed up writing, will be limited by pysam write speed
            else:
                threads = self.threads
            samtools_args = ['samtools', 'view', '-h', '-@', "%s" % threads, '-b', '/dev/stdin', '-o', self.path]
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
                sort_bam(inpath=self.path, output=self.path, sort_order=self.sort_order, threads=self.threads)
        except Exception:
            # If no reads had been written to self.path
            # an exception will be raised by is_file_coordinate_sorted.
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

    def __init__(self, path, external_bin='choose_best', sort_order='coordinate', threads=4):
        """
        Read Bam files.

        Use this class with a contexthandler.

        :param path: Path to read bam file from.
        :param external_bin: Specify `samtools` to use samtools or `None` to use pysam for writing to `path`.
        :param sort_order: Can be `coordinate` or `queryname` and will cause the output file to sorted by this strategy.
        """
        self.path = path
        self.external_bin = external_bin
        self.sort_order = sort_order
        self.threads = threads

    @property
    def args(self):
        """Figure out in sambamba or samtools are available and return correct arguments."""
        if self.external_bin:
            if self.threads > 3:
                threads = 3  # More threads won't speed up samtools reading
            else:
                threads = self.threads
            samtools_args = ['samtools', 'view', '-@', "%s" % threads, '-h', self.path]
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
            self.path = sort_bam(inpath=self.path, output=self.path, sort_order=self.sort_order, threads=self.threads)
        if self.is_bam and self.args:
            self.proc = subprocess.Popen(self.args, stdout=subprocess.PIPE, env=os.environ.copy(), close_fds=True)
            self.af = pysam.AlignmentFile(self.proc.stdout)
        else:
            self.af = pysam.AlignmentFile(self.path)
        return self.af

    def __exit__(self, type, value, traceback):
        """Provide context handler exit."""
        self.close()
