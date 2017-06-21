import copy
import gzip
import logging
import os
import shutil
import subprocess
import tempfile
import pysam
import six

if six.PY2:
    import shutilwhich  # noqa: F401

__VERSION__ = '0.3.18'
logger = logging.getLogger(__name__)
logging.basicConfig(format='%(asctime)s %(name)s %(levelname)s - %(message)s', level=logging.DEBUG)


def is_file_coordinate_sorted(path, reads_to_check=1000):
    """Determine if first 1000 reads are coordinate sorted."""
    with pysam.AlignmentFile(path) as f:
        i = 0
        current_start = 0
        current_tid = 0
        for r in f:
            if i > reads_to_check:
                return True
            if not r.is_unmapped:
                i += 1
                if r.reference_start >= current_start and r.tid >= current_tid:
                    current_start = r.reference_start
                    current_tid = r.tid
                    continue
                else:
                    logger.info("%s is not sorted by coordinate", path)
                    return False
    logger.info("%s is sorted by coordinate", path)
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
    qn = ''
    count = 0
    for r in f:
        current_query_name = r.query_name
        if current_query_name != qn:
            # Next query_name
            count += 1
            if count % chunk_size == 0 and qn:
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
    # The first read returned (even after seeking!) is always the first read in the file.
    # Calling next once resolves that.
    if not f.tell() == start:
        next(f)
        f.seek(start)
    reads = []
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
    """Return start positions that returns the first read after the current last qname."""
    f = pysam.AlignmentFile(fn)
    start = f.tell()
    last_qnames = copy.deepcopy(last_qnames)
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


def merge_bam(bam_collection, template_bam, output_path):
    """Merge a readname sorted collection of BAM files."""
    bam_collection = [bam for bam in bam_collection if os.path.exists(bam)]
    args = ['samtools', 'cat', '-h', template_bam, '-o', output_path]
    args.extend(bam_collection)
    subprocess.call(args, env=os.environ.copy())
    [shutil.rmtree(bam, ignore_errors=True) for bam in bam_collection]
    return output_path


def sort_bam(inpath, output, sort_order, threads=1):
    """Sort bam file at inpath using sort_order and write output to output."""
    fd = None
    if inpath == output:
        # We sort 'in place'
        fd, temp_out = tempfile.mkstemp()
    else:
        temp_out = output
    args = ['samtools', 'sort', '-@', "%s" % threads]
    if sort_order == 'queryname':
        args.append('-n')
    args.extend(['-o', temp_out, inpath])
    logger.info("Sorting bam file with command '%s'", " ".join(args))
    subprocess.call(args, env=os.environ.copy())
    if temp_out != output:
        shutil.move(temp_out, output)
    if fd:
        os.close(fd)
    return output


def index_bam(inpath):
    """Index BAM file at input."""
    if not os.path.exists("%s.bai" % inpath):
        pysam.index(inpath)


def split_locations_between_clusters(bamfile, self_tag='AD', other_tag='BD', distance=1000000):
    """Split a bam file into multiple chunks, where chunks have no tagged reads close to the split."""
    # TODO: modify this so we don't split small chromosomes, instead we should "merge" small chromosomes
    index_bam(bamfile)
    with pysam.AlignmentFile(bamfile) as f:
        name_length = [(chrom['SN'], chrom['LN']) for chrom in f.header['SQ']]
        chunks = []
        for (name, length) in name_length:
            chunk = []
            split_pos = [i for i in range(0, length, distance)]
            split_pos.append(length)
            start_end_pos = [list(z) for z in zip(split_pos[:-1], split_pos[1:])]
            for i, (start, end) in enumerate(start_end_pos):
                if not end == length:  # i.e this is not the last chunk
                    end = find_end(f, chrom=name, end=end, self_tag=self_tag, other_tag=other_tag)
                    start_end_pos[i + 1][0] = end  # Update the next start position with the new end
                chunk.append("%s:%s-%s" % (name, start, end))
            chunks.extend(chunk)
    return chunks


def find_end(f, chrom, end, self_tag, other_tag, padding=5000):
    """Find a position where the distance between tags is high and we can split safely."""
    min_end = end - padding
    max_end = end + padding
    current_end = min_end
    tagged_reads = [r for r in f.fetch(reference=chrom, start=min_end, end=max_end) if not r.is_unmapped and (r.has_tag(self_tag) or r.has_tag(other_tag))]
    if not tagged_reads:
        return end
    distance_to_tag = tagged_reads[0].reference_start - min_end
    pos = int((min_end + tagged_reads[0].reference_start) / 2.0)
    for r in tagged_reads:
        current_distance = r.reference_start - current_end
        if current_distance > distance_to_tag:
            distance_to_tag = current_distance
            pos = int((r.reference_start + current_end) / 2.0)
        current_end = r.reference_start
    if max_end - current_end > distance_to_tag:
        pos = int((max_end + current_end) / 2.0)
    return pos


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
        self.external_bin = None  # external_bin -- check if failing writes can be prevented by using pysam
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

    def __init__(self, path, external_bin='choose_best', sort_order=None, threads=4, region=None, index=False):
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
        self.region = region
        self.index = index
        self.setup()

    def setup(self):
        """Index and sort alignment file if necessary."""
        if self.sort_order:
            if is_file_coordinate_sorted(self.path):
                sort_order = 'coordinate'
            else:
                sort_order = 'queryname'
            if sort_order != self.sort_order:
                self.path = sort_bam(inpath=self.path, output=self.path, sort_order=self.sort_order, threads=self.threads)
        if self.region or self.index:
            index_bam(self.path)

    @property
    def args(self):
        """Figure out in sambamba or samtools are available and return correct arguments."""
        if self.external_bin:
            if self.threads > 3:
                threads = 3  # More threads won't speed up samtools reading
            else:
                threads = self.threads
            samtools_args = ['samtools', 'view', '-@', "%s" % threads, '-h', self.path]
            if self.region:
                samtools_args.append(self.region)
            return samtools_args
        return None

    @property
    def is_bam(self):
        """Return whether file is BAM (True) or SAM."""
        if not hasattr(self, '_is_bam'):
            try:
                g = gzip.GzipFile(self.path)
                self._is_bam = g.read(3) == 'BAM' or b'BAM'
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
        if self.is_bam and self.args:
            self.proc = subprocess.Popen(self.args, stdout=subprocess.PIPE, env=os.environ.copy(), close_fds=True)
            self.af = pysam.AlignmentFile(self.proc.stdout)
        else:
            self.af = pysam.AlignmentFile(self.path)
        return self.af

    def __exit__(self, type, value, traceback):
        """Provide context handler exit."""
        self.close()
