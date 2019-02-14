import copy
import logging
import os
import shutil
import tempfile
import pysam
import compare_reads

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


def get_queryname_positions(fn, chunk_size=10000, threads=0):
    """
    Get positions in order to split alignment file into equal-sized portions.

    Return a list of tuples, where the first item is the current file position for the first read,
    and the second item is the last query_name of the chunk.
    """
    with pysam.AlignmentFile(fn, threads=threads) as f:
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
            if (count + 1) % chunk_size == 0:
                last_pos = f.tell()
        seek_positions.append((start, qn))
    return seek_positions


def get_reads(fn, start, last_qname, threads=0):
    """Get reads starting at `start` and ending with last_qname."""
    with pysam.AlignmentFile(fn, threads=threads) as f:
        # The first read returned (even after seeking!) is always the first read in the file.
        # Calling next once resolves that.
        if not f.tell() == start:
            next(f)
            f.seek(start)
        reads = []
        for r in f:
            if compare_reads._compare_sort_with_queryname(r, last_qname) < 1:
                reads.append(r)
            else:
                break
    return reads


def start_positions_for_last_qnames(fn, last_qnames, threads=0):
    """Return start positions that returns the first read after the current last qname."""
    with pysam.AlignmentFile(fn, threads=threads) as f1, pysam.AlignmentFile(fn, threads=threads) as f2:
        last_qnames = copy.deepcopy(last_qnames)
        seek_positions = []
        last_seek = f1.tell()
        i = 0
        next(f1)
        current_segment = next(f1)
        next(f2)
        while last_qnames:
            current_last_qname = last_qnames.pop(0)
            try:
                while compare_reads._compare_sort_with_queryname(current_segment, current_last_qname) < 1:
                    # IOW while we haven't iterated past the expected read position
                    current_segment = next(f1)
                    next(f2)
                i += 1
                logging.info("Got chunk %i", i)
                seek_positions.append(last_seek)
                last_seek = f2.tell()
            except StopIteration:
                break
        seek_positions.append(last_seek)
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
        self.af = pysam.AlignmentFile(self.path, mode="wb", template=self.template, header=self.header, threads=self.threads)
        return self.af

    def __exit__(self, type, value, traceback):
        """Provide context handler exit."""
        self.close()


class BamAlignmentReader(object):
    """Wraps pysam.AlignmentFile with sambamba for reading if input file is a bam file."""

    def __init__(self, path, sort_order=None, threads=4, region=None, index=False):
        """
        Read Bam files.

        Use this class with a contexthandler.

        :param path: Path to read bam file from.
        :param external_bin: Specify `samtools` to use samtools or `None` to use pysam for writing to `path`.
        :param sort_order: Can be `coordinate` or `queryname` and will cause the output file to sorted by this strategy.
        """
        self.path = path
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

    def close(self):
        """Close filehandles and subprocess safely."""
        self.af.close()

    def __enter__(self):
        """Provide context handler entry."""
        self.af = pysam.AlignmentFile(self.path, threads=self.threads - 1)
        return self.af

    def __exit__(self, type, value, traceback):
        """Provide context handler exit."""
        self.close()
