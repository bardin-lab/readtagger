import pysam
from cached_property import cached_property
from .genotype import Genotype
from .cap3 import Cap3Assembly
from .tagcluster import TagCluster

import warnings


class Cluster(list):
    """A Cluster of reads."""

    left_blast_result = None
    right_blast_result = None

    def __init__(self):
        """Setup Cluster object."""
        super(Cluster, self).__init__()
        self.nref = 0
        self.id = -1

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
    def read_index(self):
        """Index of read names in cluster."""
        if not hasattr(self, '_read_index') or (hasattr(self, '_clusterlen') and len(self) != self._clusterlen):
            self._read_index = set([r.query_name for r in self])
        return self._read_index

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

    @cached_property
    def left_support(self):
        """Number of supporting reads on the left side of cluster."""
        return self.clustertag.left_sequence_count

    @cached_property
    def right_support(self):
        """Number of supporting reads on the right side of cluster."""
        return self.clustertag.right_sequence_count

    @cached_property
    def score(self):
        """Sum of all supporting reads for this cluster."""
        return self.left_support + self.right_support

    @cached_property
    def nalt(self):
        """Return number of unique read names that support an insertion."""
        return len(self.read_index)

    def _make_contigs(self):
        for contigs in self.left_contigs:
            pass
        for contigs in self.right_contigs:
            pass

    @cached_property
    def left_contigs(self):
        """Left contigs for this cluster."""
        if self.clustertag.left_sequences:
            return [contig.sequence for contig in self.clustertag.left_insert.assembly.contigs]
        else:
            return []

    @cached_property
    def right_contigs(self):
        """Right contigs for this cluster."""
        if self.clustertag.right_sequences:
            return [contig.sequence for contig in self.clustertag.right_insert.assembly.contigs]
        else:
            return []

    @cached_property
    def start(self):
        """Start coordinate for this cluster."""
        return self._start_and_end[0]

    @cached_property
    def end(self):
        """End coordinate for this cluster."""
        return self._start_and_end[1]

    @cached_property
    def _start_and_end(self):
        start = self.clustertag.five_p_breakpoint
        end = self.clustertag.three_p_breakpoint
        if start is None:
            start = end
        if end is None:
            end = start
        if start > end:
            end, start = start, end
        if start == end:
            end += 1
        return start, end

    @cached_property
    def valid_tsd(self):
        """Current cluster is a TSD."""
        return self.clustertag.tsd.is_valid

    def genotype_likelihood(self):
        r"""
        Calculate genotype likelihood for current cluster.

        P(g|D) = P(g)P(D\g)/sum(P(g)P(D|g')) where P(D|g) = Pbin(Nalt, Nalt + Nfef)
        :return:
        """
        return Genotype(nref=self.nref, nalt=self.nalt)

    def set_id(self, id):
        """Set a numeric id that identifies this cluster."""
        self.id = id

    def to_fasta(self):
        """Write contigs (or reads if no contig) out as fasta items."""
        fasta_items = []
        if self.left_contigs:
            for i, contig in enumerate(self.left_contigs):
                fasta_items.append(">cluster_%s_left_contig_%s\n%s\n" % (self.id, i, contig))
        else:
            for key, seq in self.clustertag.left_sequences.items():
                fasta_items.append(">cluster_%s_left_sequences_%s\n%s\n" % (self.id, key, seq))
        if self.right_contigs:
            for i, contig in enumerate(self.right_contigs):
                fasta_items.append(">cluster_%s_right_contig_%s\n%s\n" % (self.id, i, contig))
        else:
            for key, seq in self.clustertag.right_sequences.items():
                fasta_items.append(">cluster_%s_left_sequences_%s\n%s\n" % (self.id, key, seq))
        return fasta_items

    def serialize(self):
        """Return id, start, end and read_index for multiprocessing."""
        return (self.id, self.start, self.end, self.read_index.copy())


def non_evidence(data):
    """Count all reads that point against evidence for a transposon insertion."""
    result = {}
    input_path = data['input_path']
    chromosome = data['chromosome']
    chunk = data['chunk']
    first_chunk = chunk[0]
    last_chunk = chunk[-1]
    start = first_chunk[1]
    min = start - 500
    if min <= 0:
        min = 1
    end = last_chunk[2]
    max = end + 500
    try:
        f = pysam.AlignmentFile(input_path)
        try:
            reads = f.fetch(chromosome, min, max)
        except Exception:
            pysam.index(input_path)
            f = pysam.AlignmentFile(input_path)
            reads = f.fetch(chromosome, min, max)
        for r in reads:
            if not r.is_duplicate and r.is_proper_pair and r.mapq > 0:
                add_to_clusters(chunk, r, result)
        f.close()
    except ValueError:
        warnings.warn("Encountered ValueError on chromosome %s for start %s and end %s of chunks %s" % (chromosome, start, end, chunk))
    return result


def add_to_clusters(chunk, r, result):
    """Manage non-evidence results."""
    if r.is_supplementary:
        min_start = r.pos
        max_end = r.aend
    else:
        min_start = min([r.pos, r.mpos])
        max_end = max(r.aend, r.pos + r.isize)
    for index, start, end, read_index in chunk:
        if (min_start < start < max_end and min_start < end < max_end) and r.query_name not in read_index:
            # A read is only incompatible if it overlaps both ends
            read_index.add(r.query_name)  # We count fragments only once
            if index not in result:
                result[index] = 1
            else:
                result[index] += 1
