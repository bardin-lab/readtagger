import os
import tempfile
from cached_property import cached_property
from concurrent.futures import (
    wait,
    ThreadPoolExecutor
)

from .bam_io import BamAlignmentReader as Reader
from .bam_io import BamAlignmentWriter as Writer
from .gff_io import write_cluster
from .tagcluster import TagCluster
from .cap3 import Cap3Assembly
from .verify import discard_supplementary


class Cluster(list):
    """A Cluster of reads."""

    @cached_property
    def min(self):
        """
        Cache leftmost start of cluster.

        This assumes that the cluster is filled from left to right.
        """
        return self[0].pos

    @cached_property
    def tid(self):
        """Cache current reference id."""
        return self[0].tid

    @property
    def max(self):
        """Reference end of last read added to cluster."""
        return self[-1].reference_end

    def overlaps(self, r):
        """Determine if r overlaps the current cluster."""
        return r.pos <= self.max

    def same_chromosome(self, r):
        """Whether r is on same chromsome as cluster."""
        return self.tid == r.tid

    def read_is_compatible(self, r):
        """Determine if read overlaps cluster and is on same chromosome."""
        return self.overlaps(r) and self.same_chromosome(r)

    @property
    def hash(self):
        """Calculate a hash based on read name and read sequence for all reads in this cluster."""
        string_to_hash = "|".join(["%s%s" % (r.query_name, r.query_sequence) for r in self])
        return hash(string_to_hash)

    @property
    def clustertag(self):
        """Return clustertag for current cluster of reads."""
        if not hasattr(self, '_clustertag') or (hasattr(self, '_clusterlen') and len(self) != self._clusterlen):
            self._clusterlen = len(self)
            self._clustertag = TagCluster(self)
        return self._clustertag

    def can_join(self, other_cluster, max_distance=1500):
        """
        Join clusters that have been split or mates that are not directly connected.

        This can happen when an insertion erodes a number of nucleotides at the insertion site. Example in HUM4 at chr3R:13,373,242-13,374,019.
        The situation looks like this then:

        Reference Genome:       012345678901234567890
        left cluster read1:     >>>>>>>>XXXXX
        left cluster read2:       >>>>>>XXXXX
        right cluster read1             XXXXX<<<<<<<<
        right cluster read2:            XXXXX<<<<<

        Bunch of characteristics here:
          - left cluster shouldn't have any 3'support, right cluster no 5' support
          - clipped reads should overlap (except if a large number of nucleotides have been eroded: that could be an interesting mechanism.)
          - inferred insert should point to same TE # TODO: implement this
        """
        if hasattr(self, '_cannot_join_d'):
            # We already tried joining this cluster with another cluster,
            # so if we checked if we can try joining the exact same clusters
            # (self.hash is key in self._cannot_join and other_cluster.hash is value in self._cannot_join)
            # we know this didn't work and save ourselves the expensive assembly check
            other_hash = self._cannot_join_d.get(self.hash)
            if other_hash == other_cluster.hash:
                return False
        if hasattr(self, '_can_join_d'):
            other_hash = self._can_join_d.get(self.hash)
            if other_hash == other_cluster.hash:
                return True
        return self._can_join(other_cluster, max_distance)

    def _can_join(self, other_cluster, max_distance):
        other_clustertag = TagCluster(other_cluster)
        # First check ... are three_p and five_p of cluster overlapping?
        if not self.clustertag.tsd.three_p and not other_clustertag.tsd.five_p:
            if self.clustertag.tsd.five_p and other_clustertag.tsd.three_p:
                extended_three_p = other_clustertag.tsd.three_p - other_clustertag.tsd.three_p_clip_length
                extended_five_p = self.clustertag.tsd.five_p_clip_length + self.clustertag.tsd.five_p
                if extended_three_p <= extended_five_p:
                    self._can_join_d = {self.hash: other_cluster.hash}
                    return True
        # Next check ... can informative parts of mates be assembled into the proper insert sequence
        if self.clustertag.left_sequences and other_clustertag.left_sequences:
            # A cluster that provides support for a 5p insertion will have the reads always annotated as left sequences.
            # That's a bit confusing, since the mates are to the right of the cluster ... but that's how it is.
            if (other_clustertag.five_p_breakpoint - self.clustertag.five_p_breakpoint) < max_distance:
                # We don't want clusters to be spaced too far away. Not sure if this is really a problem in practice.
                if Cap3Assembly.sequences_contribute_to_same_contig(self.clustertag.left_sequences, other_clustertag.left_sequences):
                    self._can_join_d = {self.hash: other_cluster.hash}
                    return True
        # We know this cluster (self) cannot be joined with other_cluster, so we cache this result,
        # Since we may ask this question multiple times when joining the clusters.
        self._cannot_join_d = {self.hash: other_cluster.hash}
        return False


class ClusterFinder(object):
    """Find clusters of reads."""

    def __init__(self, input_path,
                 output_bam=None,
                 output_gff=None,
                 include_duplicates=False,
                 sample_name=None,
                 threads=1,
                 min_mapq=1,
                 max_clustersupport=200,
                 remove_supplementary_without_primary=True):
        """
        Find readclusters in input_path file.

        This class assumes that all reads in input_path are potentially interesting and that the alignment has been done for paired end reads.
        You will not want to run this on a full high coverage alignment file, since the clusters will become huge.
        The initial limits of the cluster are defined by the maximum extent of overlapping reads, and each read that is added at the 3 prime end of
        the cluster will extend the cluster.
        The join_cluster method will then join clusters that overlap through their clipped sequences and cluster that can be assembled based on their proximity
        and the fact that they support the same same insertion (and can hence contribute to the same contig if assembled).
        """
        self.input_path = input_path
        self._sample_name = sample_name
        self.output_bam = output_bam
        self.output_gff = output_gff
        self.include_duplicates = include_duplicates
        self.min_mapq = min_mapq
        self.max_clustersupport = max_clustersupport
        self.remove_supplementary_without_primary = remove_supplementary_without_primary
        self.threads = threads
        self.tp = ThreadPoolExecutor(threads)  # max theads
        self.cluster = self.find_cluster()
        self.clean_clusters()
        self.join_clusters()
        self.to_bam()
        self.to_gff()

    @cached_property
    def sample_name(self):
        """Return sample name if passed in manually, else guess sample name from input file."""
        if not self._sample_name:
            basename = os.path.basename(self.input_path)
            if '.' in basename:
                basename = basename.rsplit('.', 1)[0]
            return basename
        else:
            return self._sample_name

    def _remove_supplementary_without_primary(self):
        """Remove supplementary reads without primary alignments."""
        output_path = tempfile.mktemp(dir='.')
        discard_supplementary(input_path=self.input_path, output_path=output_path)
        self.input_path = output_path

    def find_cluster(self):
        """Find clusters by iterating over input_path and creating clusters if reads are disjointed."""
        if self.remove_supplementary_without_primary:
            self._remove_supplementary_without_primary()
        clusters = []
        with Reader(self.input_path) as reader:
            self.header = reader.header
            for r in reader:
                if not self.include_duplicates:
                    if r.is_duplicate:
                        continue
                if not r.mapping_quality > self.min_mapq:
                    continue
                if not clusters:
                    cluster = Cluster()
                    cluster.append(r)
                    clusters.append(cluster)
                    continue
                if clusters[-1].read_is_compatible(r):
                    clusters[-1].append(r)
                else:
                    cluster = Cluster()
                    cluster.append(r)
                    clusters.append(cluster)
        return clusters

    def clean_clusters(self):
        """Remove clusters that have more reads supporting an insertion than specified in self.max_clustersupport."""
        self.cluster = [c for c in self.cluster if not sum([len(c.clustertag.left_sequences), len(c.clustertag.right_sequences)]) > self.max_clustersupport]

    def join_clusters(self):
        """Iterate over self.cluster and attempt to join consecutive clusters."""
        if len(self.cluster) > 1:
            cluster_length = len(self.cluster)
            new_clusterlength = 0
            while new_clusterlength != cluster_length:
                # Iterate until the length of the cluster doesn't change anymore.
                self._cache_join()
                cluster_length = new_clusterlength
                prev_cluster = self.cluster[0]
                for cluster in self.cluster[1:]:
                    if prev_cluster.can_join(cluster):
                        prev_cluster.extend(cluster)
                        self.cluster.remove(cluster)
                    else:
                        prev_cluster = cluster
                new_clusterlength = len(self.cluster)

    def _cache_join(self):
        futures = []
        prev_cluster = self.cluster[0]
        for cluster in self.cluster[1:]:
            futures.append(self.tp.submit(prev_cluster.can_join, cluster))
            prev_cluster = cluster
        wait(futures)

    def to_bam(self):
        """Write clusters of reads and include cluster number in CD tag."""
        if self.output_bam:
            with Writer(self.output_bam, header=self.header) as writer:
                for i, cluster in enumerate(self.cluster):
                    for r in cluster:
                        r.set_tag('CD', i)
                        writer.write(r)

    def to_gff(self):
        """Write clusters as GFF file."""
        if self.output_gff:
            write_cluster(self.cluster, self.header, self.output_gff, sample=self.sample_name, threads=self.threads)
