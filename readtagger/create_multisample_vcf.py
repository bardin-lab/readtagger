from collections import deque
import pysam


from .cluster import Cluster
from .find_softclip_clusters import SoftClipCluster
from .edlib_align import multiple_sequences_overlap
from .filter_insertions import (
    sequences_match,
    SCAN_SOFT_CLIP_REGION
)
from .utils import overlap

# SEARCH_WINDOW is the maximum distance around an insertion in which to look for potentially overlapping insertions
SEARCH_WINDOW = 800
VCF_MANDATORY = list(Cluster.vcf_mandatory.keys())
VCF_ME_INFO = list(Cluster.vcf_info.keys())
VCF_ME_SAMPLE = list(Cluster.vcf_sample.keys())
VCF_SC_INFO = list(SoftClipCluster.vcf_info.keys())
VCF_SC_SAMPLE = list(SoftClipCluster.vcf_sample.keys())


class NotSingleSampleVcfException(Exception):
    """Indicate that function can't deal with multi sample VCF files."""


def window(variant_files, chrom, n=1000):
    """Create a sliding window of size n using collections.deque."""
    it = yield_sorted_records(variant_files, chrom=chrom)
    win = deque((next(it, None) for _ in range(n)), maxlen=n)
    yield win
    for e in it:
        win.append(e)
        yield win


def collect_sample_names(variant_files):
    """Return the sample names given an iterable of VariantFile objects."""
    sample_names = []
    for f in variant_files:
        sample_name = list(f.header.samples)
        if len(sample_name) != 1:
            raise NotSingleSampleVcfException('Can only process single sample VCF/BCF files, offending file is "%s"' % f.filename)
        else:
            sample_names.append(sample_name[0])
    return sample_names


class VCFMerger(object):
    """Class that will merge multiple VCF samples."""

    def __init__(self, variant_file_paths, output_path, window_size=1000, search_window=SEARCH_WINDOW):
        """Merge VCF file in variant_file_paths and write merged output to output_path."""
        self.variant_file_paths = variant_file_paths
        self.variant_files = [pysam.VariantFile(f) for f in self.variant_file_paths]
        self.sample_names = collect_sample_names(self.variant_files)
        self.header = self.setup_header()
        self.window_size = window_size
        self.search_window = search_window
        with pysam.VariantFile(output_path, 'w', header=self.header) as self.merged_variants:
            self.evaluate()

    def setup_header(self):
        """Add new samples to header."""
        header = self.variant_files[0].header.copy()
        for sample_name in self.sample_names[1:]:
            header.add_sample(sample_name)
        return header

    def evaluate(self):
        """Process all variant files."""
        def process_chunks(chunks):
            chunk = chunks.popleft()
            if chunk:
                merged_record = self.merge_items(current_record=chunk, records=chunks)
                self.fix_empty_values(merged_record)
                self.merged_variants.write(merged_record)

        for chrom in list(self.header.contigs):
            windows = window(self.variant_files, chrom=chrom, n=self.window_size)
            for chunks in windows:
                process_chunks(chunks)
            while chunks:
                process_chunks(chunks)

    def fix_empty_values(self, record):
        """Remove consecutive `.` in empty format fields."""
        # This works around a bug in pysam which seems to fill missing values to the longest samples' value
        # This sounds similar to what BCF encoders/decoders are supposed to do.
        for sample, vrs in record.samples.items():
            for key, value in vrs.items():
                if value is not None and isinstance(value, str):
                    split_value = value.split()
                    if split_value and '\x07' in split_value[0]:
                        record.samples[sample][key] = None

    def merge_items(self, current_record, records):
        """Merge all compatible items into current record and drop compatible item from records."""
        # Need to create a new record using the new header containing all samples
        current_name = next(iter(current_record.samples))
        current_record = self.copy_record(current_record=current_record)
        for j, other_record in enumerate(records):
            if other_record:
                if current_record.start > other_record.stop + self.search_window:
                    break
                else:
                    if self.can_merge(current_record, other_record, current_name):
                        current_record = self.merge_record(current_record, other_record)
                        records[j] = None
        return current_record

    def can_merge(self, current_record, other_record, current_name):
        """Determine if current_record and other_record can be merged."""
        current_svtype = current_record.info['SVTYPE']
        other_svtype = other_record.info['SVTYPE']
        current_start = current_record.start
        other_start = other_record.start
        if current_svtype.startswith('SOFTCLIP') and current_svtype == other_svtype:
            if abs(current_start - other_start) < SCAN_SOFT_CLIP_REGION:
                # We assume there's only 1 consensus
                current_sequence = current_record.samples[current_name]['CLIP_CONSENSUS']
                other_sequence = [vrs['CLIP_CONSENSUS'] for vrs in other_record.samples.values()][0]
                compare = '5p_clip' if current_svtype == 'SOFTCLIP:5P' else '3p_clip'
                if sequences_match(current_sequence, other_sequence, compare=compare):
                    return True
        elif current_svtype == other_svtype == 'INS:ME':
            if overlap(start1=current_start, end1=current_record.stop, start2=other_start, end2=other_record.stop):
                current_mename = self.get_format_values(record=current_record, key='MENAME')
                other_mename = self.get_format_values(record=other_record, key='MENAME')
                if current_mename & other_mename:
                    return True
                else:
                    for key in ('MEASSEMBLY5', 'MEASSEMBLY3'):
                        current_assembled_sequences = self.get_format_values(record=current_record, key=key)
                        if current_assembled_sequences:
                            other_assembled_sequences = self.get_format_values(record=other_record, key=key)
                            if other_assembled_sequences and multiple_sequences_overlap(current_assembled_sequences,
                                                                                        other_assembled_sequences,
                                                                                        check_revcom=False):
                                return True
                # TODO: if MENAME doesn't match try to rescue via annotated softclips
        return False

    @staticmethod
    def get_format_values(record, key):
        """Return all FORMAT values for all sample of a VCF record."""
        values = set()
        for sample in record.samples.values():
            value = sample[key]
            if value:
                if isinstance(value, tuple):
                    for v in value:
                        if v and len(v) > 10:
                            values.add(v)
                else:
                    values.add(value)
        return values

    def merge_record(self, current_record, other_record):
        """Update current_record with format keys from other_record."""
        # Pre-requisite is that current record has been copied for the new headers
        other_name = next(iter(other_record.samples.keys()))
        # TODO: Merge some of the INFO keys (MATE_IDs is probably useful?)
        # We copy the other records' format attributes
        for k, v in other_record.samples[other_name].items():
            current_record.samples[other_name][k] = v
        current_record.stop = min((other_record.stop, current_record.stop))
        current_record.start = max((other_record.start, current_record.start))
        return current_record

    def copy_record(self, current_record):
        """Create a new record using an updated header and copy attributes."""
        new_record = self.header.new_record()
        # Now copy the "static" attributes for this record
        for key in VCF_MANDATORY:
            setattr(new_record, key, getattr(current_record, key))
        if current_record.info['SVTYPE'] == 'INS:ME':
            VCF_INFO = VCF_ME_INFO
        else:
            VCF_INFO = VCF_SC_INFO
        for key in VCF_INFO:
            new_record.info[key] = current_record.info[key]
        current_name = next(iter((current_record.samples.keys())))
        for k, v in current_record.samples[current_name].items():
            new_record.samples[current_name][k] = v
        return new_record


def yield_sorted_records(variant_files, chrom='2L'):
    """Yield from a list of variant files records in sorted order."""
    current_stack = {}
    iterators = [f.fetch(contig=chrom) for f in variant_files]
    for i, vf in enumerate(iterators):
        try:
            record = next(vf)
            current_stack[i] = record
        except StopIteration:
            pass
    while current_stack:
        (key, record) = min(current_stack.items(), key=lambda x: x[1].start)
        yield record
        try:
            current_stack[key] = next(iterators[key])
        except StopIteration:
            del current_stack[key]
