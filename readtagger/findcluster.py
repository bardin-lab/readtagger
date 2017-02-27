from cached_property import cached_property

from .bam_io import BamAlignmentReader as Reader
from .tagcluster import TagCluster


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

    def can_join(self, other_cluster):
        """
        Join clusters that have been split.

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
        self_clustertag = TagCluster(self)
        other_clustertag = TagCluster(other_cluster)
        # TODO: definitely need to check if clipped sequences overlap, or introduce maximum distance between clusters
        if not self_clustertag.tsd.three_p and not other_clustertag.tsd.five_p:
            if self_clustertag.tsd.five_p and other_clustertag.tsd.three_p:
                extended_three_p = other_clustertag.tsd.three_p - other_clustertag.tsd.three_p_clip_length
                extended_five_p = self_clustertag.tsd.five_p_clip_length + self_clustertag.tsd.five_p
                if extended_three_p <= extended_five_p:
                    return True
                else:
                    return False
        return False


class ClusterFinder(object):
    """Find clusters of reads."""

    def __init__(self, input_path):
        """
        Find readclusters in input_path file.

        This class assumes that all reads in input_path are potentially interesting and that the alignment has been done for paired end reads.
        You will not want to run this on a full high coverage alignment file, since the clusters will become huge.
        The initial limits of the cluster are defined by the maximum matepair distance, and each read that is added at the 3 prime end of
        the cluster will extend the cluster
        """
        self.input_path = input_path
        self.cluster = self.find_cluster()
        self.join_clusters()

    def find_cluster(self):
        """Find clusters by iterating over input_path and creating clusters if reads are disjointed."""
        clusters = []
        with Reader(self.input_path) as reader:
            r = next(reader)
            cluster = Cluster()
            cluster.append(r)
            clusters.append(cluster)
            for r in reader:
                if clusters[-1].read_is_compatible(r):
                    clusters[-1].append(r)
                else:
                    cluster = Cluster()
                    cluster.append(r)
                    clusters.append(cluster)
        return clusters

    def join_clusters(self):
        """Iterate over self.cluster and attempt to join consecutive clusters."""
        if len(self.cluster) > 1:
            prev_cluster = self.cluster[0]
            for cluster in self.cluster[1:]:
                if prev_cluster.can_join(cluster):
                    prev_cluster.extend(cluster)
                    self.cluster.remove(cluster)
                else:
                    prev_cluster = cluster
