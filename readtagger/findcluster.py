import os
import shutil
import tempfile
from cached_property import cached_property
from concurrent.futures import (
    wait,
    ThreadPoolExecutor
)

from .bam_io import BamAlignmentReader as Reader
from .bam_io import BamAlignmentWriter as Writer
from .cluster import Cluster
from .gff_io import write_cluster
from .verify import discard_supplementary


class ClusterFinder(object):
    """Find clusters of reads."""

    def __init__(self,
                 input_path,
                 output_bam=None,
                 output_gff=None,
                 include_duplicates=False,
                 sample_name=None,
                 threads=1,
                 min_mapq=1,
                 max_clustersupport=200,
                 remove_supplementary_without_primary=False):
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
        self._tempdir = tempfile.mkdtemp(dir='.')
        self.remove_supplementary_without_primary = remove_supplementary_without_primary
        self.threads = threads
        self.tp = ThreadPoolExecutor(threads)  # max theads
        self.cluster = self.find_cluster()
        self.clean_clusters()
        self.join_clusters()
        self.to_bam()
        self.to_gff()
        shutil.rmtree(self._tempdir, ignore_errors=True)

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
        output_path = os.path.join(self._tempdir, 'clean.bam')
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
