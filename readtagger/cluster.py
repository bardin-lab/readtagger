from collections import (
    defaultdict,
    OrderedDict
)
from hashlib import md5
from itertools import (
    chain,
    groupby,
    permutations
)
from .edlib_align import (
    multiple_sequences_overlap,
    sequences_overlap
)
from .cap3 import Cap3Assembly
from .genotype import Genotype
from .instance_lru import instance_method_lru_cache
from .tagcluster import TagCluster
from .utils import overlap

MIN_LONG_READ = 200
# Minimum read length to consider a read coming from longread tech
MAX_VALID_ISIZE = 10000
# Should cover a lot of deletions, while sufficiently small to not include reads
# that aligned to a TE somewhere else on the same chromosome
MIN_VALID_ISIZE_FOR_NON_PROPER_PAIR = 700
# If a read is not a proper pair we still consider it if is not a proper pair because
# of an isize that is too big, which can happen with deletions.
MAX_COLLECT_EVIDENCE = 10000
# maximum amount of evidence to consider


class BaseCluster(list):
    """Common attributes for clusters of reads."""

    def __init__(self):
        """Initialize BaseCluster instance."""
        super(BaseCluster, self).__init__()
        self.nref = 0
        self.evidence_against = set()
        self.evidence_for_five_p = set()
        self.evidence_for_three_p = set()
        self.feature_args = []
        self.exclude = False
        self.ref = 'N'
        self.sequence = -1

    def __hash__(self):
        """Delegate to self.hash for hash specific to this cluster."""
        return self.hash

    def __eq__(self, other):
        """Define equality as mathcing hashes."""
        return self.hash == other.hash

    def __ne__(self, other):
        """Define not equal as not equal."""
        return not self == other

    @property
    def reference_name(self):
        """Return current reference name."""
        if self:
            return next(iter(self)).reference_name

    @property
    def tid(self):
        """Cache current reference id."""
        return next(iter(self)).tid

    @property
    def read_index(self):
        """Index of read names in cluster."""
        return set(r.query_name for r in chain(self, self.evidence_for_five_p, self.evidence_for_three_p) if not r.has_tag('AC'))  # AC is assembled contig

    @property
    def hash(self):
        """Calculate a hash based on read name and read sequence for all reads in this cluster."""
        return hash(tuple((id(r) for r in self)))

    @property
    def max_mapq(self):
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

    @property
    def genotype(self):
        """Return most likely genotype for this cluster."""
        return self.genotype_likelihoods.genotype

    @property
    def genotype_likelihoods(self):
        r"""
        Calculate genotype likelihood for current cluster.

        P(g|D) = P(g)P(D\g)/sum(P(g)P(D|g')) where P(D|g) = Pbin(Nalt, Nalt + Nfef)
        :return:
        """
        return Genotype(nref=self.nref, nalt=self.nalt)


class Cluster(BaseCluster):
    """A Cluster of reads."""

    exportable = ['source', 'score', 'total_left_count', 'left_mate_count',
                  'total_right_count', 'right_mate_count', 'nref', 'proper_pair_only',
                  'max_mapq', 'genotype', 'genotype_likelihoods', 'valid_TSD',
                  'left_inserts', 'right_inserts', 'insert_reference_name', 'softclip_clusters']
    vcf_mandatory = OrderedDict([
        ('alts', 'alts',),
        ('pos', 'pos',),
        ('stop', 'stop',),
        ('chrom', 'reference_name',),
        ('ref', 'ref',),
        ('id', 'vcf_id',),
    ])
    vcf_info = OrderedDict([
        ('SVTYPE', 'svtype',),
        ('SVLEN', 'insert_len',),
        ('MQ', 'max_mapq',),
        ('EVENT', 'vcf_id',),
        ('MATEID', 'softclip_clusters',),
        ('VALID_TSD', 'valid_tsd'),
    ])
    vcf_sample = OrderedDict([
        ('GT', 'vcf_genotype',),
        ('GL', 'vcf_genotype_likelikoods',),
        ('AD', ['nref', 'nalt'],),
        ('DP', 'depth',),
        ('SU', 'nalt',),
        ('SU5', 'total_left_count',),
        ('SU3', 'total_right_count',),
        ('PE', 'total_mate_count'),
        ('PE5', 'left_mate_count'),
        ('PE3', 'right_mate_count'),
        ('SR', 'total_split_count',),
        ('SR5', 'left_split_count',),
        ('SR3', 'right_split_count',),
        ('MSP', 'evidence_spanning_insertion'),
        ('MENAME', 'insert_reference_name',),
        ('MESTART', 'insert_start',),
        ('MEEND', 'insert_end',),
        ('MEASSEMBLY5', 'left_inserts',),
        ('MEASSEMBLY3', 'right_inserts',),
    ])
    source = "findcluster"

    def __init__(self, shm_dir, max_proper_size=0):
        """Initialize Cluster instance."""
        super(Cluster, self).__init__()
        self.insert_reference_name = None
        self.max_proper_size = max_proper_size
        self.shm_dir = shm_dir
        self._cannot_join_d = {}
        self.abnormal = False
        self.softclip_clusters = []

    @property
    def proper_pair_only(self):
        """Indicate if cluster is composed of reads that are exclusively aligned as proper pairs."""
        for r in self:
            if not r.is_proper_pair:
                return False
        return True

    @property
    def evidence_spanning_insertion(self):
        """Return the number of mate pair fragments that span 5' and 3' of an insertion."""
        return len(set(self.clustertag.left_sequences.keys()) & set(self.clustertag.right_sequences.keys()))

    @property
    def pos(self):
        """Return the 1-based start position of this cluster."""
        return self.start + 1

    @property
    def vcf_genotype(self):
        """Return the genotype in VCF format."""
        gt = {'homozygous': (1, 1),
              'heterozygous': (0, 1),
              'reference': (0, 0)}
        return gt[self.genotype]

    @property
    def vcf_genotype_likelikoods(self):
        """Return genitype likelihood for VCF output."""
        return tuple(self.genotype_likelihoods)

    @property
    def vcf_id(self):
        """Return a valid VCF ID."""
        return str(self.id)

    @property
    def alts(self):
        """Return alts for VCF output."""
        return ('<INS:ME>',)

    @property
    def svtype(self):
        """Return the SVTYPE for VCF output."""
        return 'INS:ME'

    @property
    def id(self):
        """Return a unique id for this cluster."""
        unique_string = "{reference_name},{sequence},{start}".format(reference_name=self.reference_name,
                                                                     sequence=self.sequence,
                                                                     start=self.start)
        return "INS_%s" % md5(unique_string.encode()).hexdigest()

    @property
    def stop(self):
        """Return the stop for VCF output."""
        return self.end + 1

    @property
    def depth(self):
        """Return an estimate of the total depth."""
        # TODO: this is just a first pass, the correct thing would be to include uninformative reads
        return self.nref + self.nalt

    @property
    def total_split_count(self):
        """Return number of split reads supporting this insertion."""
        return self.left_split_count + self.right_split_count

    @property
    def left_split_count(self):
        """Return a count of all evidence for the 5prime of an insertion."""
        return self.total_left_count - self.left_mate_count

    @property
    def right_split_count(self):
        """Return a count of all evidence for the 5prime of an insertion."""
        return self.total_right_count - self.right_mate_count

    @property
    def insert_start(self):
        """Return approximate start coordinate of insert."""
        if self.feature_args:
            return self.feature_args[0].sbjct_start
        return None

    @property
    def insert_end(self):
        """Return approximate end coordinate of insert."""
        if self.feature_args:
            return self.feature_args[0].sbjct_end - self.feature_args[0].sbjct_start
        return None

    @property
    def insert_len(self):
        """Return best guess of reference insert length."""
        if self.feature_args:
            return self.feature_args[0].sbjct_end - self.feature_args[0].sbjct_start
        return None

    @property
    def type(self):
        """Return the type of insert."""
        return self.insert_reference_name or "TE"

    @property
    def min(self):
        """
        Cache leftmost start of cluster.

        This assumes that the cluster is filled from left to right.
        """
        return min((r.reference_start for r in self))

    @property
    def max(self):
        """Return reference end of last read added to cluster."""
        return max((r.reference_end for r in self))

    def overlaps(self, r, strict=False):
        """Determine if r overlaps the current cluster."""
        if strict:
            # We use strict mode when refining cluster members, as they are not necessarily added
            # from low reference_start to high_reference_start
            start1 = min((_.reference_start for _ in self))
            end1 = max((_.reference_end for _ in self))
            return overlap(start1, end1, r.reference_start, r.reference_end)
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
        """Check that mate orientation within cluster is consistent."""
        if read.has_tag('BD'):
            # We got a mate, this means it cannot extend either 3' or 5' of the breakpoint
            if read.is_reverse:
                if read.reference_start < self.start:
                    return False
            else:
                if read.reference_end > self.end:
                    return False
        if self.clustertag.tsd.is_valid and read.has_tag('AD'):
            if read.reference_start > self.end_corrected + 10:
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
                if self.clustertag.tsd.is_valid or other_cluster.clustertag.tsd.is_valid:
                    check = Cluster(shm_dir=self.shm_dir)
                    check.extend(self + other_cluster)
                    if not check.clustertag.tsd.is_valid:
                        continue
                self.extend(other_cluster)
                self.softclip_clusters.extend(other_cluster.softclip_clusters)
                all_clusters.remove(other_cluster)

    def reachable(self, all_clusters):
        """Find all cluster that are closeby."""
        idx = all_clusters.index(self)
        current_end = self.end_corrected or self.end
        for i, cluster in enumerate(all_clusters[idx + 1:]):
            other_start = cluster.start_corrected or cluster.start
            if other_start and current_end and other_start <= current_end:
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
    def _mark_clusters_incompatible(*clusters):
        # We know these clusters have been split on purpose, don't try to merge them back together!
        for (cluster_a, cluster_b) in permutations(clusters, r=2):
            cluster_a._cannot_join_d[cluster_b.hash] = cluster_a.hash

    def check_cluster_consistency(self):
        """
        Check that clusters are internally consistent.

        If we have left or right mates without AD tags the breakpoint cannot be within the region covered by the mates.
        """
        ALIGNMENT_ERROR = 5  # Occasionally ends should be clipped, but due to mirohomologies this doesn't always happen
        if self.abnormal:
            return [self]
        three_p_reads_to_to_discard = set()
        five_p_reads_to_to_discard = set()
        for read in self.right_mate_support.values():
            if self.clustertag.three_p_breakpoint and read.reference_start < self.clustertag.three_p_breakpoint - ALIGNMENT_ERROR:
                for support_read in self.clustertag.tsd.three_p_reads:
                    three_p_reads_to_to_discard.add(support_read)
        for read in self.left_mate_support.values():
            if self.clustertag.five_p_breakpoint and read.reference_end > self.clustertag.five_p_breakpoint + ALIGNMENT_ERROR:
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
            # we know this didn't work and save ourselves the expensive assembly check,
            # but we can check if the insert reference name is compatible now
            return self.check_overlap_compatible(other_cluster, max_distance)
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
                    return True
        # Next check ... can informative parts of mates be assembled into the proper insert sequence
        if self.clustertag.left_sequences and other_clustertag.left_sequences:
            # A cluster that provides support for a 5p insertion will have the reads always annotated as left sequences.
            # That's a bit confusing, since the mates are to the right of the cluster ... but that's how it is.
            if (other_clustertag.five_p_breakpoint - self.clustertag.five_p_breakpoint) < max_distance:
                # We don't want clusters to be spaced too far away. Not sure if this is really a problem in practice.
                if multiple_sequences_overlap(self.clustertag.left_sequences.values(), other_clustertag.left_sequences.values()):
                    return True
        # TODO: Refactor this to a common function for left-left and right-right assembly
        if self.clustertag.right_sequences and other_clustertag.right_sequences:
            # A cluster that provides support for a 5p insertion will have the reads always annotated as left sequences.
            # That's a bit confusing, since the mates are to the right of the cluster ... but that's how it is.
            if (other_clustertag.three_p_breakpoint - self.clustertag.three_p_breakpoint) < max_distance:
                # We don't want clusters to be spaced too far away. Not sure if this is really a problem in practice.
                if multiple_sequences_overlap(self.clustertag.right_sequences.values(), other_clustertag.right_sequences.values()):
                    return True
        if self.check_overlap_compatible(other_cluster, max_distance):
            return True
        # We know this cluster (self) cannot be joined with other_cluster, so we cache this result,
        # Since we may ask this question multiple times when joining the clusters.
        self._mark_clusters_incompatible(self, other_cluster)
        return False

    def check_overlap_compatible(self, other_cluster, max_distance):
        """
        Check that overlap between clusters is compatible with joining.

        This means the orientation switches make sense (i.e F self and R other is ok,
        but not R self and F other), and that the putative insert reference is compatible,
        or that insert tags match (this is from the tagging step).
        """
        self_switches = self.orientation_switches
        other_switches = other_cluster.orientation_switches
        single_orientation = len(self_switches) == len(other_switches) == 1
        tolerance = 0 if self.valid_tsd else max_distance
        if abs(other_cluster.min - self.max) < max_distance and overlap(start1=self.start,
                                                                        end1=self.end,
                                                                        start2=other_cluster.start,
                                                                        end2=other_cluster.end,
                                                                        tolerance=tolerance):
            if single_orientation:
                # This check doesn't take into account the read length, as the isize
                # is determined as read-length 1 +  inner distance  + readlength 2
                min_read_read_length = min([r.reference_length for r in other_cluster if r.reference_start == other_cluster.min])
                max_read_read_length = min([r.reference_length for r in self if r.reference_end == self.max])
                if abs(other_cluster.min - self.max) < (max_distance - min_read_read_length - max_read_read_length):
                    # Merge a cluster that is split by an insertion and not connected via split reads
                    if self_switches[0][0] == 'F':
                        if set(self.insert_reference_tags()) & set(other_cluster.insert_reference_tags()):
                            return True
            elif len(other_switches) == 1 and other_switches[0][0] == 'R':
                # This is not a very specific check, so we require the most common reference of both clusters to check
                if set(self.insert_reference_tags()) & set(other_cluster.insert_reference_tags()):
                    return True

    def insert_reference_tags(self):
        """Return insert reference tags."""
        insert_ref_tags = []
        for r in self:
            if self.insert_reference_name and self.insert_reference_name != 'TE':
                insert_ref_tags.append(self.insert_reference_name)
                # This should be authoritative IMO
                continue
            if r.has_tag('AR'):
                insert_ref_tags.append(r.get_tag('AR'))
            if r.has_tag('BR'):
                insert_ref_tags.append(r.get_tag('BR'))
        return insert_ref_tags

    @property
    def total_left_count(self):
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
    def total_mate_count(self):
        """Return the number of all mates that support this cluster."""
        return self.left_mate_count + self.right_mate_count

    @property
    def total_right_count(self):
        """Return number of supporting reads on the right side of cluster."""
        return self.clustertag.right_sequence_count + len(self.evidence_for_three_p)

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
    def left_inserts(self):
        """Get a list of assembled inserts, prefixed with index."""
        return self.left_contigs

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
    def right_inserts(self):
        """Get a list of assembled inserts, prefixed with index."""
        return self.right_contigs

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
        if start is None and end is None:
            start = self.start
            end = self.end
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

    def to_fasta(self):
        """Write contigs (or reads if no contig) out as fasta items."""
        fasta_items = []
        if self.left_contigs:
            for i, contig in enumerate(self.left_contigs):
                fasta_items.append(">%s|lcontigs|%s\n%s\n" % (self.id, i, contig))
        else:
            for key, seq in self.clustertag.left_sequences.items():
                fasta_items.append(">%s|lsequences|%s\n%s\n" % (self.id, key, seq))
        if self.right_contigs:
            for i, contig in enumerate(self.right_contigs):
                fasta_items.append(">%s|rcontigs|%s\n%s\n" % (self.id, i, contig))
        else:
            for key, seq in self.clustertag.right_sequences.items():
                fasta_items.append(">%s|lsequences|%s\n%s\n" % (self.id, key, seq))
        return fasta_items

    def serialize(self):
        """Return id, start, end and read_index for multiprocessing."""
        bp_sequences = defaultdict(set)
        five_p_breakpoint = self.clustertag.five_p_breakpoint
        left_sequence = self.clustertag.left_breakpoint_sequence
        if five_p_breakpoint and left_sequence:
            bp_sequences[five_p_breakpoint] = {left_sequence}
        three_p_breakpoint = self.clustertag.three_p_breakpoint
        right_sequence = self.clustertag.right_breakpoint_sequence
        if three_p_breakpoint and right_sequence:
            bp_sequences[three_p_breakpoint].add(right_sequence)
        single_breakpoint = None
        if five_p_breakpoint and three_p_breakpoint:
            single_breakpoint = False
        elif five_p_breakpoint:
            single_breakpoint = five_p_breakpoint
        elif three_p_breakpoint:
            single_breakpoint = three_p_breakpoint
        return self.start, self.end, bp_sequences, single_breakpoint


def collect_evidence(cluster, alignment_file):
    """Count all reads that point against evidence for a transposon insertion."""
    chromosome = cluster.reference_name
    start = cluster.start
    end = cluster.end
    min_start = start - 500
    if min_start < 0:
        # Avoid pysam error for negative start coordinates
        min_start = 0
    max_end = end + 500
    start, end, bp_sequence, single_breakpoint = cluster.serialize()
    reads = alignment_file.fetch(chromosome, min_start, max_end)
    for i, r in enumerate(reads):
        if i <= MAX_COLLECT_EVIDENCE:
            if not r.is_duplicate \
                and r.mapq > 0 \
                and (r.is_proper_pair or
                     r.next_reference_name == r.reference_name == chromosome and
                     MIN_VALID_ISIZE_FOR_NON_PROPER_PAIR > abs(r.isize) < MAX_VALID_ISIZE):
                add_to_clusters(cluster, r, start, end, bp_sequence, single_breakpoint)
    cluster.evidence_against = {r for r in cluster.evidence_against if r.query_name not in cluster.read_index}
    cluster.nref = len(set(r.query_name for r in cluster.evidence_against))


def add_to_clusters(cluster, r, start, end, bp_sequence, single_breakpoint):
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

    query_name = r.query_name
    if query_name not in cluster.read_index:
        if bp_sequence:
            if {reference_start, reference_end} & set(bp_sequence):
                evidence = evidence_for(read=r, breakpoint_sequences=bp_sequence)
                if evidence:
                    if evidence == 'five_p':
                        cluster.evidence_for_five_p.add(r)
                    else:
                        cluster.evidence_for_three_p.add(r)
                    return
        if single_breakpoint:
            # We only know where one of the breakpoints is, so we ask if any reads overlap that breakpoint
            # If there is actually a TSD but we didn't detect it we slightly bias the count in favor of counting against
            # the insertion.
            # TODO: this doesn't work if we actually don't know the breakpoint, i.e when a cluster is defined
            # only by mate pairs. In that instance it might be more accurate to sample the coverage at the breakpoint
            # boundaries, and assume that  nalt / coverage estimates the AF.
            if (min_start + 1 < single_breakpoint < max_end - 1):
                cluster.evidence_against.add(r)
        elif end - start < 50:
            if min_start + 1 < start < max_end - 1 and min_start + 1 < end < max_end - 1:
                # A read is only incompatible if it overlaps both ends
                # We require the overlap to be more than 1 nucleotide (by adding 1 to min_start and subtracting 1 from max_end)
                # to avoid dealing with reads with a single mismatch at the start/end,
                # which wouldn't be soft-clipped. This shouldn't introduce any bias since we also can't assign these
                # reads to an insertion, so we simple ignore them.
                cluster.evidence_against.add(r)
        else:
            # We were not able to narrow down the insertion breakpoints.
            # We can estimate the insertion frequency by looking at how many reads overlap
            # start and end of the insertion. This isn't very precise, but insertions without
            # exact start/end are probably low in frequency anyways.
            if (min_start + 1 < start < max_end - 1) or (min_start + 1 < end < max_end - 1):
                cluster.evidence_against.add(r)


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
