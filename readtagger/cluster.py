from collections import defaultdict
from itertools import (
    chain,
    groupby,
    permutations
)

import pysam
from .edlib_align import (
    multiple_sequences_overlap,
    sequences_overlap
)
from .cap3 import Cap3Assembly
from .genotype import Genotype
from .instance_lru import instance_method_lru_cache
from .tagcluster import TagCluster

MIN_LONG_READ = 200
# Minimum read length to consider a read coming from longread tech
MAX_VALID_ISIZE = 10000
# Should cover a lot of deletions, while sufficiently small to not include reads
# that aligned to a TE somewhere else on the same chromosome
MIN_VALID_ISIZE_FOR_NON_PROPER_PAIR = 700
# If a read is not a proper pair we still consider it if is not a proper pair because
# of an isize that is too big, which can happen with deletions.


class Cluster(list):
    """A Cluster of reads."""

    left_blast_result = None
    right_blast_result = None

    def __init__(self, shm_dir, max_proper_size=0):
        """Initialize Cluster object."""
        super(Cluster, self).__init__()
        self.nref = 0
        self.id = -1
        self.feature_args = None
        self.max_proper_size = max_proper_size
        self.reference_name = None
        self.shm_dir = shm_dir
        self.evidence_against = set()
        self.evidence_for_five_p = set()
        self.evidence_for_three_p = set()
        self._can_join_d = {}
        self._cannot_join_d = {}
        self.abnormal = False

    def __hash__(self):
        """Delegate to self.hash for hash specific to this cluster."""
        return self.hash

    def __eq__(self, other):
        """Define equality as mathcing hashes."""
        return self.hash == other.hash

    def __ne__(self, other):
        """Define not equal as not equal."""
        return self != other

    @property
    def min(self):
        """
        Cache leftmost start of cluster.

        This assumes that the cluster is filled from left to right.
        """
        return min((r.reference_start for r in self))

    @property
    def tid(self):
        """Cache current reference id."""
        return self[0].tid

    @property
    def max(self):
        """Return reference end of last read added to cluster."""
        return max((r.reference_end for r in self))

    @staticmethod
    def _strict_overlap(start1, end1, start2, end2):
        """Check that range (start1, end1) overlaps with (start2, end2)."""
        # Taken from https://nedbatchelder.com/blog/201310/range_overlap_in_two_compares.html
        return end1 >= start2 and end2 >= start1

    def overlaps(self, r, strict=False):
        """Determine if r overlaps the current cluster."""
        if strict:
            # We use strict mode when refining cluster members, as they are not necessarily added
            # from low reference_start to high_reference_start
            start1 = min((_.reference_start for _ in self))
            end1 = max((_.reference_end for _ in self))
            return self._strict_overlap(start1, end1, r.reference_start, r.reference_end)
        return r.reference_start <= self.max

    def same_chromosome(self, r):
        """Whether r is on same chromsome as cluster."""
        return self.tid == r.tid

    def read_is_compatible(self, r, strict=False):
        """
        Ensure read is compatible with  cluster.

        Determines if read
          - overlaps cluster
          - is on same chromosome
          - read orientation is consistent with this cluster.
        """
        return self.overlaps(r, strict=strict) and self.same_chromosome(r) and self.read_consistent_with_clusters(r)

    def read_consistent_with_clusters(self, read):
        """Check that mater orientation within cluster is consistent."""
        if read.has_tag('BD'):
            # We got a mate, this means it cannot extend either 3' or 5' of the breakpoint
            if read.is_reverse:
                if read.reference_start < self.start:
                    return False
            else:
                if read.reference_end > self.end:
                    return False
        return True

    @property
    def orientation_switches(self):
        """Return all orientation switches in this cluster."""
        return [next(group[1]) for group in groupby(self.orientation_vector, key=lambda x: x[0])]

    def refine_members(self, assembly_realigner):
        """Try to recover more reads that support a specific insertion."""
        if assembly_realigner and not self.abnormal:
            informative_reads = assembly_realigner.collect_reads(self)
            for read in informative_reads:
                if self.read_is_compatible(read, strict=True):
                    self.append(read)

    def join_adjacent(self, all_clusters):
        """Join clusters that can be joined."""
        for other_cluster in self.reachable(all_clusters=all_clusters):
            if not self.abnormal and self.can_join(other_cluster, max_distance=self.max_proper_size):
                self.extend(other_cluster)
                all_clusters.remove(other_cluster)

    def reachable(self, all_clusters):
        """Find all cluster that are closeby."""
        idx = all_clusters.index(self)
        current_end = self.end_corrected
        for i, cluster in enumerate(all_clusters[idx + 1:]):
            if cluster.start_corrected and current_end and cluster.start_corrected <= current_end:
                yield cluster
            elif i == 0:
                yield cluster
            else:
                break

    def split_cluster_at_polarity_switch(self):
        """
        Split cluster if the direction of the mates switches.

        This is quite a rough estimate, should probaby do some more checks to ensure a single
        stray read in the wrong orientation does not split a cluster.
        """
        if self.abnormal:
            return (None, None)
        switches = self.orientation_switches
        putative_break = None
        # Make sure
        if len(switches) > 2 and switches[0][0] == 'F':
            putative_break = switches[2][1]
        if len(switches) > 2 and switches[0][0] == 'R':
            # Clusters shouldn't really start with reverse BD reads
            putative_break = switches[1][1]
        if putative_break:
            cluster_a, cluster_b = self._make_new_clusters()
            cluster_a.extend(self[:putative_break])
            cluster_b.extend(self[putative_break:])
            cluster_a, cluster_b = self.assign_reads_to_split(cluster_a, cluster_b)
            self._mark_clusters_incompatible(cluster_a, cluster_b)
            return (cluster_a, cluster_b)
        return (None, None)

    def _make_new_clusters(self, count=2):
        """Return a set of new clusters."""
        return [Cluster(shm_dir=self.shm_dir, max_proper_size=self.max_proper_size) for _ in range(count)]

    @staticmethod
    def _mark_clusters_compatible(*clusters):
        # We know these clusters have been split on purpose, don't try to merge them back together!
        for (cluster_a, cluster_b) in permutations(clusters, r=2):
            cluster_a._can_join_d[cluster_b.hash] = cluster_a.hash

    @staticmethod
    def _mark_clusters_incompatible(*clusters):
        # We know these clusters have been split on purpose, don't try to merge them back together!
        for (cluster_a, cluster_b) in permutations(clusters, r=2):
            cluster_a._cannot_join_d[cluster_b.hash] = cluster_a.hash

    def check_cluster_consistency(self):
        """
        Check that clusters are internally consistent.

        If we have left or right mates without AD tags the breakpoint cannot be within the region covered by the mates.
        """
        if self.abnormal:
            return [self]
        three_p_reads_to_to_discard = set()
        five_p_reads_to_to_discard = set()
        for read in self.right_mate_support.values():
            if self.clustertag.three_p_breakpoint and read.reference_start < self.clustertag.three_p_breakpoint:
                for support_read in self.clustertag.tsd.three_p_reads:
                    three_p_reads_to_to_discard.add(support_read)
        for read in self.left_mate_support.values():
            if self.clustertag.five_p_breakpoint and read.reference_end > self.clustertag.five_p_breakpoint:
                for support_read in self.clustertag.tsd.three_p_reads:
                    five_p_reads_to_to_discard.add(support_read)
        new_clusters = [self]
        if three_p_reads_to_to_discard:
            new_three_p_cluster = Cluster(shm_dir=self.shm_dir, max_proper_size=self.max_proper_size)
            for read in three_p_reads_to_to_discard:
                new_three_p_cluster.append(read)
                self.remove(read)
            new_clusters.append(new_three_p_cluster)
        if five_p_reads_to_to_discard:
            new_five_p_cluster = Cluster(shm_dir=self.shm_dir, max_proper_size=self.max_proper_size)
            for read in five_p_reads_to_to_discard:
                new_five_p_cluster.append(read)
                if read not in three_p_reads_to_to_discard:
                    # Can't remove read twice ...
                    # I suppose this can happen with very small TE inserts where we can have mates on both sides
                    self.remove(read)
            new_clusters.append(new_five_p_cluster)
        self._mark_clusters_incompatible(*new_clusters)
        # TODO: may want to remove any mates that correspond to the removed reads
        return new_clusters

    def assign_reads_to_split(self, putative_cluster_a, putative_cluster_b):
        """
        Check reads to split using an assembly strategy.

        Before really splitting a cluster, assemble all contributing reads and make sure there isn't a single stray read that causes
        a perfectly valid cluster to be split. If only the read at the polarity switch does not contribute to the contigs keep all other reads.
        """
        all_reads = {}
        all_reads.update(self.clustertag.left_sequences)
        all_reads.update(self.clustertag.right_sequences)
        assembly = Cap3Assembly(all_reads)
        contigs = assembly.contigs
        contig_reads = []
        cluster_a_contigs = set()
        # Establish a list of contigs and their readnames,
        # And classify whether contigs belong to cluster a or cluster b.
        for index, contig in enumerate(contigs):
            contig_reads.append(set())
            for read in contig.reads:
                readname = read.rd.name.rstrip('.1').rstrip('.2')
                if readname in putative_cluster_a.read_index:
                    cluster_a_contigs.add(index)
                contig_reads[index].add(readname)
        # collapse cluster contig reads into a or b
        cluster_a_contig_reads = set()
        for index in cluster_a_contigs:
            cluster_a_contig_reads.update(contig_reads[index])
        contig_sequences = {contig.sequence for i, contig in enumerate(contigs) if i in cluster_a_contigs}
        if contig_sequences:
            # Not all reads are assigned to contigs,
            # so we use bwa to check if a read belongs to a contig
            if len(putative_cluster_a.orientation_switches) > 1:
                switch_reads = putative_cluster_a[putative_cluster_a.orientation_switches[1][1]:]
                # Are the reads that caused the orientation switch actually contributing to cluster a?
                # Because if they don't, they should probably be ignored and treated as a separate insertion.
                sequences_to_align = [r.get_tag('MS') for r in switch_reads if r.has_tag('MS')]
                if sequences_to_align and not multiple_sequences_overlap(queries=sequences_to_align, targets=contig_sequences):
                    putative_cluster_a = Cluster(shm_dir=self.shm_dir, max_proper_size=self.max_proper_size)
                    putative_cluster_a.extend(switch_reads)
                    putative_cluster_b = Cluster(shm_dir=self.shm_dir, max_proper_size=self.max_proper_size)
                    putative_cluster_b.extend(r for r in self if r.query_name not in putative_cluster_a.read_index)
                    return putative_cluster_a, putative_cluster_b
            reads_to_remove = set()
            for read in putative_cluster_b:
                if read.query_name in cluster_a_contig_reads:
                    putative_cluster_a.append(read)
                    reads_to_remove.add(read)
                elif read.has_tag('MS') and read.has_tag('BD'):
                    ms = read.get_tag('MS')
                    for i, contig in enumerate(contig_sequences):
                        if i in cluster_a_contigs and sequences_overlap(query=ms, targets=contig_sequences):
                            putative_cluster_a.append(read)
                            reads_to_remove.add(read)
            for read in reads_to_remove:
                putative_cluster_b.remove(read)
        return putative_cluster_a, putative_cluster_b

    @property
    def orientation_vector(self):
        """Return orientation of all reads with 'BD' tag."""
        return self._get_orientation_vector()

    @instance_method_lru_cache(maxsize=10000)
    def _get_orientation_vector(self):
        return [('R', i) if r.is_reverse else ('F', i) for i, r in enumerate(self) if not r.has_tag('AD')]

    @property
    def read_index(self):
        """Index of read names in cluster."""
        return set(r.query_name for r in chain(self, self.evidence_for_five_p, self.evidence_for_three_p) if not r.has_tag('AC'))  # AC is assembled contig

    @property
    def hash(self):
        """Calculate a hash based on read name and read sequence for all reads in this cluster."""
        return hash("".join(("%s" % id(r) for r in self)))

    @property
    def clustertag(self):
        """Return clustertag for current cluster of reads."""
        return self._get_clustertag()

    @instance_method_lru_cache(maxsize=10000)
    def _get_clustertag(self):
        return TagCluster(self, shm_dir=self.shm_dir)

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
        if self.abnormal or other_cluster.abnormal:
            return False
        self_hash = self.hash
        other_hash = other_cluster.hash
        if self._cannot_join_d.get(other_hash, None) == self_hash:
            # We already tried joining this cluster with another cluster,
            # so if we checked if we can try joining the exact same clusters
            # (self.hash is key in self._cannot_join and other_cluster.hash is value in self._cannot_join)
            # we know this didn't work and save ourselves the expensive assembly check
            return False
        elif self._can_join_d.get(other_hash, None) == self_hash:
            return True
        return self._can_join(other_cluster, max_distance)

    def _can_join(self, other_cluster, max_distance):
        other_clustertag = other_cluster.clustertag
        # TODO: there should be no additional polarity switch if clusters are to be joined?
        # Second check ... are three_p and five_p of cluster overlapping?
        if not self.clustertag.tsd.three_p and not other_clustertag.tsd.five_p:
            if self.clustertag.tsd.five_p and other_clustertag.tsd.three_p:
                extended_three_p = other_clustertag.tsd.three_p - other_clustertag.tsd.three_p_clip_length
                extended_five_p = self.clustertag.tsd.five_p_clip_length + self.clustertag.tsd.five_p
                if extended_three_p <= extended_five_p:
                    self._mark_clusters_compatible(self, other_cluster)
                    return True
        # Next check ... can informative parts of mates be assembled into the proper insert sequence
        if self.clustertag.left_sequences and other_clustertag.left_sequences:
            # A cluster that provides support for a 5p insertion will have the reads always annotated as left sequences.
            # That's a bit confusing, since the mates are to the right of the cluster ... but that's how it is.
            if (other_clustertag.five_p_breakpoint - self.clustertag.five_p_breakpoint) < max_distance:
                # We don't want clusters to be spaced too far away. Not sure if this is really a problem in practice.
                if multiple_sequences_overlap(self.clustertag.left_sequences.values(), other_clustertag.left_sequences.values()):
                    self._mark_clusters_compatible(self, other_cluster)
                    return True
        # TODO: Refactor this to a common function for left-left and right-right assembly
        if self.clustertag.right_sequences and other_clustertag.right_sequences:
            # A cluster that provides support for a 5p insertion will have the reads always annotated as left sequences.
            # That's a bit confusing, since the mates are to the right of the cluster ... but that's how it is.
            if (other_clustertag.three_p_breakpoint - self.clustertag.three_p_breakpoint) < max_distance:
                # We don't want clusters to be spaced too far away. Not sure if this is really a problem in practice.
                if multiple_sequences_overlap(self.clustertag.right_sequences.values(), other_clustertag.right_sequences.values()):
                    self._mark_clusters_compatible(self, other_cluster)
                    return True
        self_switches = self.orientation_switches
        other_switches = other_cluster.orientation_switches
        if len(self_switches) == 1 and len(other_switches) == 1 and self_switches != other_switches:
            if abs(other_cluster.min - self.max) < max_distance and self._strict_overlap(start1=self.start,
                                                                                         end1=self.end,
                                                                                         start2=other_cluster.start,
                                                                                         end2=other_cluster.end):
                # This check doesn't take into account the read length, as the isize
                # is determined as read-length 1 +  inner distance  + readlength 2
                min_read_read_length = min([r.reference_length for r in other_cluster if r.reference_start == other_cluster.min])
                max_read_read_length = min([r.reference_length for r in self if r.reference_end == self.max])
                if abs(other_cluster.min - self.max) < (max_distance - min_read_read_length - max_read_read_length):
                    # Merge a cluster that is split by an insertion and not connected via split reads
                    if self_switches[0][0] == 'F' and other_switches[0][0] == 'R':
                        self._mark_clusters_compatible(self, other_cluster)
                        return True
        # We know this cluster (self) cannot be joined with other_cluster, so we cache this result,
        # Since we may ask this question multiple times when joining the clusters.
        self._mark_clusters_incompatible(self, other_cluster)
        return False

    @property
    def total_left_support(self):
        """Return number of supporting reads on the left side of cluster."""
        return self.clustertag.left_sequence_count + len(self.evidence_for_five_p)

    @property
    def left_mate_support(self):
        """
        Return mates on the left that support an insertion.

        This is excluding split reads.
        """
        return {r.query_name: r for r in self if r.has_tag('BD') and not r.has_tag('AD') and
                ("%s.1" % r.query_name in self.clustertag.left_sequences or "%s.2" % r.query_name in self.clustertag.left_sequences)}

    @property
    def left_mate_count(self):
        """Return number of mates on the left that support an insertion."""
        return len(self.left_mate_support)

    @property
    def right_mate_support(self):
        """
        Return mates on the right that support an insertion.

        This is excluding split reads.
        """
        return {r.query_name: r for r in self if r.has_tag('BD') and not r.has_tag('AD') and
                ("%s.1" % r.query_name in self.clustertag.right_sequences or "%s.2" % r.query_name in self.clustertag.right_sequences)}

    @property
    def right_mate_count(self):
        """Return mates on the right that support an insertion."""
        return len(self.right_mate_support)

    @property
    def total_right_support(self):
        """Return number of supporting reads on the right side of cluster."""
        return self.clustertag.right_sequence_count + len(self.evidence_for_three_p)

    @property
    def maximum_mapq(self):
        """Return the highest MAPQ observed for this cluster."""
        return max((r.mapq for r in self))

    @property
    def score(self):
        """Return sum of all supporting reads for this cluster."""
        return self.nalt

    @property
    def nalt(self):
        """Return number of unique read names that support an insertion."""
        # The read index contains all read names that contribute to this cluster (proper cluster sequences, but also split reads
        # picked up with `evidence_for`)
        return len(self.read_index)

    def _make_contigs(self):
        # We just touch the contigs in threading mode,
        # so the actual assembly is being triggered.
        for contigs in self.left_contigs:
            pass
        for contigs in self.right_contigs:
            pass

    @property
    def left_contigs(self):
        """Left contigs for this cluster."""
        return self._get_left_contigs()

    @instance_method_lru_cache(maxsize=10000)
    def _get_left_contigs(self):
        if not self.abnormal and self.clustertag.left_sequences:
            return [contig.sequence for contig in self.clustertag.left_insert.contigs]
        else:
            return []

    @property
    def right_contigs(self):
        """Right contigs for this cluster."""
        return self._get_right_contigs()

    @instance_method_lru_cache(maxsize=10000)
    def _get_right_contigs(self):
        if not self.abnormal and self.clustertag.right_sequences:
            return [contig.sequence for contig in self.clustertag.right_insert.contigs]
        else:
            return []

    @property
    def start(self):
        """Start coordinate for this cluster."""
        return self._start_and_end[0]

    @property
    def end(self):
        """End coordinate for this cluster."""
        return self._start_and_end[1]

    @property
    def start_corrected(self):
        """
        Extend start site.

        May be extended by maximum insert size length minus the most distant 5' mate if no exact TSD exists.
        """
        if not self.clustertag.tsd.three_p_support and len(self.orientation_switches) < 2:
            # No TSD and no orientation switch, we don't know where the cluster starts
            three_p_reads = self.right_mate_support.values()
            if three_p_reads:
                start_position = min([r.reference_start for r in three_p_reads])
                end_position = max([r.reference_end for r in three_p_reads])
                if end_position - start_position < self.max_proper_size:
                    return end_position - self.max_proper_size
        return self.clustertag.five_p_breakpoint

    @property
    def end_corrected(self):
        """
        Extend end site.

        May be extended by maximum insert size length plus the most distant 3' mate end if no exact TSD exists.
        """
        if not self.clustertag.tsd.five_p_support and len(self.orientation_switches) < 2:
            # No TSD and no orientation switch, we don't know where the cluster ends
            three_p_reads = self.left_mate_support.values()
            if three_p_reads:
                start_position = min([r.reference_start for r in three_p_reads])
                end_position = max([r.reference_end for r in three_p_reads])
                if end_position - start_position < self.max_proper_size:
                    return start_position + self.max_proper_size
        return self.clustertag.three_p_breakpoint

    @property
    def _start_and_end(self):
        return self._get_start_and_end()

    @instance_method_lru_cache(maxsize=10000)
    def _get_start_and_end(self):
        start = self.start_corrected
        end = self.end_corrected
        if start is None:
            start = end
        if end is None:
            end = start
        if start < 0:
            start = 0
        if start > end:
            end, start = start, end
        if start == end:
            end += 1
        return start, end

    @property
    def valid_tsd(self):
        """Return whether current cluster has a valid TSD."""
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
                fasta_items.append(">cluster_%s_left_contigs_%s\n%s\n" % (self.id, i, contig))
        else:
            for key, seq in self.clustertag.left_sequences.items():
                fasta_items.append(">cluster_%s_left_sequences_%s\n%s\n" % (self.id, key, seq))
        if self.right_contigs:
            for i, contig in enumerate(self.right_contigs):
                fasta_items.append(">cluster_%s_right_contigs_%s\n%s\n" % (self.id, i, contig))
        else:
            for key, seq in self.clustertag.right_sequences.items():
                fasta_items.append(">cluster_%s_left_sequences_%s\n%s\n" % (self.id, key, seq))
        return fasta_items

    def serialize(self):
        """Return id, start, end and read_index for multiprocessing."""
        bp_sequences = {}
        five_p_breakpoint = self.clustertag.five_p_breakpoint
        left_sequence = self.clustertag.left_breakpoint_sequence
        if five_p_breakpoint and left_sequence:
            bp_sequences[five_p_breakpoint] = {left_sequence}
        three_p_breakpoint = self.clustertag.three_p_breakpoint
        right_sequence = self.clustertag.right_breakpoint_sequence
        if three_p_breakpoint and right_sequence:
            if three_p_breakpoint not in bp_sequences:
                bp_sequences[three_p_breakpoint] = {right_sequence}
            else:
                bp_sequences[three_p_breakpoint].add(right_sequence)
        if five_p_breakpoint and three_p_breakpoint:
            single_breakpoint = False
        elif five_p_breakpoint:
            single_breakpoint = five_p_breakpoint
        elif three_p_breakpoint:
            single_breakpoint = three_p_breakpoint
        else:
            single_breakpoint = None  # Shouldn't happen IRL I think
        return (self.id, self.start, self.end, self.read_index.copy(), bp_sequences, single_breakpoint)


def non_evidence(data):
    """Count all reads that point against evidence for a transposon insertion."""
    result = {'against': {},
              'for': {}}
    input_path = data['input_path']
    chromosome = data['chromosome']
    chunk = data['chunk']
    # Chunk is a small list of clusters
    start = min((c[1] for c in chunk))
    end = max((c[2] for c in chunk))
    min_start = start - 500
    if min_start < 0:
        # Avoid pysam error for negative start coordinates
        min_start = 0
    max_end = end + 500
    with pysam.AlignmentFile(input_path) as f:
        try:
            reads = f.fetch(chromosome, min_start, max_end)
        except Exception:
            pysam.index(input_path)
            f = pysam.AlignmentFile(input_path)
            reads = f.fetch(chromosome, min_start, max_end)
        for r in reads:
            if not r.is_duplicate \
                and r.mapq > 0 \
                and (r.is_proper_pair or
                     r.next_reference_name == r.reference_name == chromosome and
                     MIN_VALID_ISIZE_FOR_NON_PROPER_PAIR > abs(r.isize) < MAX_VALID_ISIZE):
                add_to_clusters(chunk, r, result)
    for index in result['against']:
        for_reads = set(result['for'][index])
        against_reads = set(result['against'][index])
        result['against'][index] = {k: result['against'][index][k] for k in against_reads - for_reads}
    return result


def add_to_clusters(chunk, r, result):
    """
    Count reads overlapping a cluster.

    If a read r overlaps a cluster region,
    but does not show evidence for an insertion it will be counted.
    Once a read name has been seen it will not be counted again.
    """
    reference_start = r.reference_start
    reference_end = r.reference_end
    if r.is_supplementary or r.alen > MIN_LONG_READ:  # supplementary or long read
        min_start = reference_start
        max_end = reference_end
    else:
        min_start = min([reference_start, r.next_reference_start])
        max_end = max(reference_end, r.reference_start + r.isize)
    # TODO: the for loop below is a little wasteful, since we call it for each read,
    # but a single read is unlikely to provide evidence for more than 2 or 3 clusters.
    # We could jump ahhead in the loop somehow.
    for index, start, end, supporting_read_index, bp_sequence, single_breakpoint in chunk:
        if index not in result['against']:
            result['against'][index] = defaultdict(list)
        if index not in result['for']:
            result['for'][index] = dict()
        query_name = r.query_name
        if query_name not in supporting_read_index and query_name not in result['for'][index]:
            if bp_sequence:
                if {reference_start, reference_end} & set(bp_sequence):
                    evidence = evidence_for(read=r, breakpoint_sequences=bp_sequence)
                    if evidence:
                        result['for'][index][query_name] = (r, evidence)
                        continue
            if single_breakpoint:
                # We only know where one of the breakpoints is, so we ask if any reads overlap that breakpoint
                # If there is actually a TSD but we didn't detect it we slightly bias the count in favor of counting against
                # the insertion.
                # TODO: this doesn't work if we actually don't know the breakpoint, i.e when a cluster is defined
                # only by mate pairs. In that instance it might be more accurate to sample the coverage at the breakpoint
                # boundaries, and assume that  nalt / coverage estimates the AF.
                if (min_start + 1 < single_breakpoint < max_end - 1):
                    result['against'][index][query_name].append(r)
            elif end - start < 50:
                if min_start + 1 < start < max_end - 1 and min_start + 1 < end < max_end - 1:
                    # A read is only incompatible if it overlaps both ends
                    # We require the overlap to be more than 1 nucleotide (by adding 1 to min_start and subtracting 1 from max_end)
                    # to avoid dealing with reads with a single mismatch at the start/end,
                    # which wouldn't be soft-clipped. This shouldn't introduce any bias since we also can't assign these
                    # reads to an insertion, so we simple ignore them.
                    result['against'][index][query_name].append(r)
            else:
                # We were not able to narrow down the insertion breakpoints.
                # We can estimate the insertion frequency by looking at how many reads overlap
                # start and end of the insertion. This isn't very precise, but insertions without
                # exact start/end are probably low in frequency anyways.
                if (min_start + 1 < start < max_end - 1) or (min_start + 1 < end < max_end - 1):
                    result['against'][index][query_name].append(r)


def evidence_for(read, breakpoint_sequences):
    """
    Check if the clipped sequence of a read supports an insertion.

    `read` is a pysam AlignedSegment, breakpoint_sequences is a dictionary,
    where keys is the breakpoint and breakpoint_sequences is the value.
    """
    # TODO: if we don't have any breakpoint sequences this will (perhaps falsely) return False
    # Probably better to underestimate this and not come up with a fancy solution
    # TODO: allow matching IUPAC letters
    bp_sequences = breakpoint_sequences.get(read.reference_end)
    if bp_sequences:
        soft_clipped_sequence = read.seq[read.qend:]
        if soft_clipped_sequence:
            if any(s.startswith(soft_clipped_sequence[:4]) for s in bp_sequences):
                return 'five_p'
    bp_sequences = breakpoint_sequences.get(read.reference_start)
    if bp_sequences:
        soft_clipped_sequence = read.seq[:read.qstart]
        if soft_clipped_sequence:
            if any(s.endswith(soft_clipped_sequence[-4:]) for s in bp_sequences):
                return 'three_p'
    return False
