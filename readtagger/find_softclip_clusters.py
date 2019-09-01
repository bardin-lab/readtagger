import logging
from hashlib import md5

from .bam_io import BamAlignmentReader as Reader
from .cluster import BaseCluster
from .cluster_base import (
    SampleNameMixin,
    ToGffMixin
)
from .dumb_consensus import dumb_consensus
from .tag_softclip import get_softclipped_portion

logger = logging.getLogger(__name__)


class SoftClipCluster(BaseCluster):
    """A cluster that groups reads with the same soft clipping position."""

    exportable = ['source', 'consensus', 'score', 'max_mapq', 'ID']
    vcf_mandatory = {
        'alts': 'alts',
        'pos': 'pos',
        'chrom': 'reference_name',
        'ref': 'ref',
        'id': 'id',
    }
    vcf_info = {
        'SVTYPE': 'vcf_type',
        'MQ': 'max_mapq',
        'EVENT': 'mate_id',
        'MATEID': 'mate_id'
    }
    vcf_sample = {
        'AD': ['nref', 'nalt'],
        'SU': 'nalt',
        'CLIP_CONSENSUS': 'consensus'
    }
    source = 'find_softclip'

    def __init__(self, clip_position, clip_type, sample_name):
        """Collect all reads with same `clip_position` and `clip_type`."""
        super(SoftClipCluster, self).__init__()
        self._sample_name = sample_name
        self.clip_position = clip_position
        self.start = self.clip_position
        self.end = self.clip_position
        self.clip_type = clip_type
        self.type = clip_type
        self.clipped_sequences = []
        self.mate_id = None
        self.nref = 0

    def set_id(self, id):
        """Set a numeric id that identifies this cluster."""
        self.id = id

    @property
    def pos(self):
        """Return the 1-based start of the softclip position."""
        return self.start + 1

    @property
    def alts(self):
        """Return VCF formatted ALT allele."""
        return ('<%s>' % self.vcf_type,)

    @property
    def vcf_type(self):
        """Return the VCF compatible type of this cluster."""
        if self.clip_type == '3p_clip':
            return 'SOFTCLIP:3P'
        return 'SOFTCLIP:5P'

    def append(self, read, seq=None):
        """Append read to `self`."""
        super(SoftClipCluster, self).append(read)
        if seq:
            self.clipped_sequences.append(seq)

    def read_is_compatible(self, clip_position, clip_type):
        """Determine if `clip_position` and `clip_type` are compatible with this cluster."""
        return clip_position == self.clip_position and clip_type == self.clip_type

    @property
    def consensus(self):
        """Return a consensus sequence for the clipped sequences im this cluster."""
        return dumb_consensus(string_list=self.clipped_sequences, left_align=self.clip_type == '3p_clip')

    def reachable(self, all_clusters, idx):
        """Find all cluster that are at the same position and have the same clip_type."""
        clip_position = self.clip_position
        clip_type = self.clip_type
        for i, cluster in enumerate(all_clusters[idx + 1:]):
            if cluster.clip_position == clip_position:
                if cluster.clip_type == clip_type:
                    yield cluster
            else:
                break


class SoftClipClusterFinder(SampleNameMixin, ToGffMixin):
    """Find clusters of softclipped reads."""

    def __init__(self, input_path=None, region=None, min_mapq=4, sample_name=None, output_gff=None, threads=1):
        """Find and report clusters of softclipped reads."""
        self._sample_name = sample_name
        self.threads = threads
        self.input_path = input_path
        self.output_gff = output_gff
        self.region = region
        self.min_mapq = min_mapq
        self.header = None
        self.clusters = []
        if input_path:
            self.find_clusters()
            self.merge_clusters()
            self.to_gff(output_gff)

    def find_clusters(self):
        """Find clusters by iterating over input_path and creating clusters if reads are softclipped."""
        logger.info("Finding clusters of softclipped reads (%s)" % self.region or 0)
        with Reader(self.input_path, index=True, sort_order='coordinate') as reader:
            self.header = reader.header
            for r in reader.fetch(region=self.region):
                if not r.is_duplicate and r.mapping_quality >= self.min_mapq:
                    self.add_read(r=r)
        logger.info("Found %d clusters (%s)", len(self.clusters), self.region or 0)

    def add_read(self, r):
        """Add a clipped read."""
        softclipped_portions = get_softclipped_portion(read=r, min_clip_length=2)
        for read, start, end in softclipped_portions:
            if start == 0:
                clip_position = r.reference_start
                clip_type = '5p_clip'
            elif end == r.query_length:
                clip_position = r.reference_end
                clip_type = '3p_clip'
            # That should cover all relevant possibilities
            seq = r.query_sequence[start:end]
            if not self.clusters:
                # This is the first cluster we found
                cluster = SoftClipCluster(clip_position=clip_position, clip_type=clip_type, sample_name=self.sample_name)
                self.clusters.append(cluster)
            else:
                cluster = self.clusters[-1]
            if not cluster.read_is_compatible(clip_position, clip_type=clip_type):
                cluster = SoftClipCluster(clip_position=clip_position, clip_type=clip_type, sample_name=self.sample_name)
                self.clusters.append(cluster)
            cluster.append(read=r, seq=seq)

    def merge_clusters(self):
        """Merge clusters with same cluster_type and same `clip_type`."""
        self.clusters.sort(key=lambda x: x.clip_position)
        for i, c in enumerate(self.clusters):
            c.set_id(i)
        new_list = []
        deleted_items = 0
        for i, cluster in enumerate(self.clusters):
            for reachable in cluster.reachable(all_clusters=self.clusters, idx=i):
                cluster.extend(reachable)
                del self.clusters[reachable.id - deleted_items]
                deleted_items += 1
            new_list.append(cluster)
        self.clusters = new_list
        for i, cluster in enumerate(self.clusters):
            unique_string = "{reference_name},{sequence},{start}".format(reference_name=cluster.reference_name,
                                                                         sequence=i,
                                                                         start=cluster.start)
            cluster.set_id("SOFTCLIP_%s_%s" % (cluster._sample_name, md5(unique_string.encode()).hexdigest()))
        logger.info("Found %s clusters after merging clusters", len(self.clusters))
