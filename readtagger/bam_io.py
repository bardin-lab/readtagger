import copy
import gzip
import logging
import os
import shutil
import subprocess
import tempfile
import pysam

logger = logging.getLogger(__name__)


def is_file_coordinate_sorted(path, reads_to_check=1000):
    """Determine if first 1000 reads are coordinate sorted."""
    with pysam.AlignmentFile(path) as f:
        i = 0
        current_start = 0
        current_tid = 0
        try:
            if f.header['HD']['SO'] == 'coordinate':
                # We trust the header, seems most sane
                logger.info("%s is sorted by coordinate", path)
                return True
        except Exception:
            pass
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


def get_mean_read_length(path, reads_to_check=1000):
    """Get mean read length for the first reads in source."""
    read_length = 0
    i = 0  # In case no reads are in source, avoids potential unbounded variable error
    with pysam.AlignmentFile(path) as source:
        for i, r in enumerate(source):
            if i >= reads_to_check:
                break
            read_length += r.query_length
    mean_read_length = int(read_length / i) if i > 0 else 0
    logger.info("Mean read length is '%s'", mean_read_length)
    return mean_read_length


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
        if r.query_name == last_qname or qname_cmp_func(r.query_name, last_qname) == 1:
            # We either reached the last qname, or we are already beyond,
            # for instance if the last_qname isn't in the file. This works
            # because the files are qname sorted.
            last_qname = r.query_name
            try:
                while r.query_name == last_qname:
                    reads.append(r)
                    r = next(f)
            except StopIteration:
                break
            return reads
        else:
            reads.append(r)
    return reads


# Copied from https://github.com/10XGenomics/supernova/blob/master/tenkit/lib/python/tenkit/bam.py,
# itself based on sort_bam.c in htslib
def qname_cmp_func(qname1, qname2):
    """Return 1 if qname1 > qname2, 0 is qname1 == gname2 and -1 of qname1 < qname2."""
    def update_char(qname, index):
        """Return current char and index."""
        index += 1
        if index < len(qname):
            c = qname[index]
        else:
            c = ''
        return c, index

    def get_ord(c):
        """Return ord if c or 0."""
        if c:
            return ord(c)
        return 0

    index_i, index_j = 0, 0
    c1 = qname1[0] if qname1 else None
    c2 = qname2[0] if qname2 else None
    while c1 and c2:
        if c1.isdigit() and c2.isdigit():
            while c1 == '0':
                c1, index_i = update_char(qname1, index_i)
            while c2 == '0':
                c2, index_j = update_char(qname2, index_j)
            while c1.isdigit() and c2.isdigit() and c1 == c2:
                c1, index_i = update_char(qname1, index_i)
                c2, index_j = update_char(qname2, index_j)
            if c1.isdigit() and c2.isdigit():
                index_k, index_l = index_i, index_j
                while c1.isdigit() and c2.isdigit():
                    c1, index_k = update_char(qname1, index_k)
                    c2, index_l = update_char(qname2, index_l)
                return 1 if c1.isdigit() else (-1 if c2.isdigit() else get_ord(qname1[index_i]) - get_ord(qname2[index_j]))
            elif c1.isdigit():
                return 1
            elif c2.isdigit():
                return -1
            elif index_i != index_j:
                return 1 if index_i < index_j else -1
        else:
            if c1 != c2:
                return get_ord(c1) - get_ord(c2)
            c1, index_i = update_char(qname1, index_i)
            c2, index_j = update_char(qname2, index_j)
    return 1 if c1 else (-1 if c2 else 0)


def start_positions_for_last_qnames(fn, last_qnames):
    """Return start positions that returns the first read after the current last qname."""
    f = pysam.AlignmentFile(fn)
    start = f.tell()
    last_qnames = copy.deepcopy(last_qnames)
    current_last_qname = last_qnames.pop(0)
    seek_positions = [start]
    for r in f:
        if r.query_name == current_last_qname or qname_cmp_func(r.query_name, current_last_qname) > 0:
            current_last_qname = r.query_name
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
            if last_qnames:
                current_last_qname = last_qnames.pop(0)
            else:
                # We've reached the end of last_qnames
                return seek_positions
    return seek_positions


def merge_bam(bam_collection, output_path, template_bam=None, sort_order=None, threads=1):
    """Merge a readname sorted collection of BAM files."""
    bam_collection = [bam for bam in bam_collection if os.path.exists(bam)]
    if not template_bam:
        template_bam = bam_collection[0]
    pysam.cat('-h', template_bam, '-o', output_path, *bam_collection)
    for bam in bam_collection:
        try:
            os.remove(bam)
        except OSError:
            pass
    if sort_order:
        sort_bam(inpath=output_path, output=output_path, sort_order=sort_order, threads=threads)
    return output_path


def sort_bam(inpath, output, sort_order, threads=1, cram=False, reference_fasta=None):
    """
    Sort bam file at inpath using sort_order and write output to `output`.

    Specify `cram=True` and supply `reference_fasta` to write output as CRAM file.
    """
    fd = None
    if inpath == output:
        # We sort 'in place'
        fd, temp_out = tempfile.mkstemp(dir=os.path.dirname(output))
    else:
        temp_out = output
    args = ["-@%s" % threads]
    if sort_order == 'queryname':
        args.append('-n')
    if cram:
        args.extend(['--output-fmt=cram,no_ref'])
    args.extend(['-o', temp_out, inpath])
    pysam.sort(*args)
    if temp_out != output:
        shutil.move(temp_out, output)
    if fd:
        os.close(fd)
    return output


def index_bam(inpath):
    """Index BAM file at input."""
    if not os.path.exists("%s.bai" % inpath):
        pysam.index(inpath)


def split_locations_between_clusters(bamfile, self_tag='AD', other_tag='BD', distance=1000000, region=None):
    """Split a bam file into multiple chunks, where chunks have no tagged reads close to the split."""
    # TODO: modify this so we don't split small chromosomes, instead we should "merge" small chromosomes
    index_bam(bamfile)
    limit_chrom = None
    limit_start = None
    limit_end = None
    if region:
        region = region.split(':')
        limit_chrom = region[0]
        limit_start, limit_end = [int(i) for i in region[1].split('-')]
    with pysam.AlignmentFile(bamfile) as f:
        name_length = [(chrom['SN'], chrom['LN']) for chrom in f.header['SQ']]
        chunks = []
        for (name, length) in name_length:
            if limit_chrom and limit_chrom != name:
                continue
            chunk = []
            split_pos = [i for i in range(1, length, distance)]
            split_pos.append(length)
            start_end_pos = [list(z) for z in zip(split_pos[:-1], split_pos[1:])]
            for i, (start, end) in enumerate(start_end_pos):
                if limit_end:
                    if end > limit_end or start < limit_start:
                        continue
                if not end == length:  # i.e this is not the last chunk
                    end = find_end(f, chrom=name, end=end, self_tag=self_tag, other_tag=other_tag)
                    start_end_pos[i + 1][0] = end  # Update the next start position with the new end
                # pysam doesn't like to start coordinates at 0, so we just add +1
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
    """Wrap pysam.AlignmentFile with samtools for multithreaded compressed writing."""

    def __init__(self, path, template=None, header=None, threads=1, sort_order='coordinate'):
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
        self.threads = threads
        self.sort_order = sort_order

    def close(self):
        """
        Close filehandles and subprocess safely.

        If necessary will sort the the file.
        """
        self.af.close()
        if os.path.exists(self.path):
            sort_order = 'coordinate' if is_file_coordinate_sorted(self.path) else 'queryname'
            if sort_order != self.sort_order:
                sort_bam(inpath=self.path, output=self.path, sort_order=self.sort_order, threads=self.threads)

    def __enter__(self):
        """Provide context handler entry."""
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
            sort_order = 'coordinate' if is_file_coordinate_sorted(self.path) else 'queryname'
            if sort_order != self.sort_order:
                self.path = sort_bam(inpath=self.path, output=self.path, sort_order=self.sort_order, threads=self.threads)
        if self.region or self.index:
            index_bam(self.path)

    @property
    def args(self):
        """Figure out if samtools is available and return correct arguments."""
        if self.external_bin:
            # More threads won't speed up samtools reading
            threads = 3 if self.threads > 3 else self.threads
            samtools_args = ['samtools', 'view', "-@%s" % threads, '-h', self.path]
            if self.region:
                samtools_args.append(self.region)
            return samtools_args

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
