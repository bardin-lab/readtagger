import logging
from itertools import groupby

import pysam
from cached_property import cached_property
from .genotype import Genotype
from .bwa import SimpleAligner
from .cap3 import Cap3Assembly
from .tagcluster import TagCluster


class Cluster(list):
    """A Cluster of reads."""

    left_blast_result = None
    right_blast_result = None

    def __init__(self, shm_dir, max_proper_size=0):
        """Setup Cluster object."""
        super(Cluster, self).__init__()
        self.nref = 0
        self.id = -1
        self.feature_args = None
        self.max_proper_size = max_proper_size
        self.reference_name = None
        self.shm_dir = shm_dir

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
    def orientation_switches(self):
        """Return all orientation switches in this cluster."""
        return [next(group[1]) for group in groupby(self.orientation_vector, key=lambda x: x[0])]

    def refine_members(self, assembly_realigner):
        """Try to recover more reads that support a specific insertion."""
        if assembly_realigner:
            informative_reads = assembly_realigner.collect_reads(self)
            if informative_reads:
                self.extend(informative_reads)

    def split_cluster_at_polarity_switch(self):
        """
        Split cluster if the direction of the mates switches.

        This is quite a rough estimate, I should probaby do some more checks to ensure a single
        stray read in the wrong orientation does not split a cluster.
        """
        switches = self.orientation_switches
        putative_break = None
        # Make sure
        if len(switches) > 2 and switches[0][0] == 'F':
            putative_break = switches[2][1]
        if len(switches) > 2 and switches[0][0] == 'R':
            # Clusters shouldn't really start with reverse BD reads
            putative_break = switches[1][1]
        if putative_break:
            cluster_a = Cluster(shm_dir=self.shm_dir, max_proper_size=self.max_proper_size)
            cluster_b = Cluster(shm_dir=self.shm_dir, max_proper_size=self.max_proper_size)
            cluster_a.extend(self[:putative_break])
            cluster_b.extend(self[putative_break:])
            return self.assign_reads_to_split(cluster_a, cluster_b)
        return None, None

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
        contigs = assembly.assembly.contigs
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
        contig_sequences = [contig.sequence for i, contig in enumerate(contigs) if i in cluster_a_contigs]
        if contig_sequences:
            # Not all reads are assigned to contigs,
            # so we use bwa to check if a read belongs to a contig
            simple_aligner = SimpleAligner(reference_sequences=contig_sequences, tmp_dir=self.shm_dir)
            if len(putative_cluster_a.orientation_switches) > 1:
                switch_reads = putative_cluster_a[putative_cluster_a.orientation_switches[1][1]:]
                # Are the reads that caused the orientation switch actually contributing to cluster a?
                # Because if they don't, they should probably be ignored and treated as a separate insertion.
                sequences_to_align = [r.get_tag('MS') for r in switch_reads if r.has_tag('MS')]
                if sequences_to_align and not simple_aligner.align(sequences_to_align):
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
                        if i in cluster_a_contigs and simple_aligner.align(ms):
                            putative_cluster_a.append(read)
                            reads_to_remove.add(read)
            for read in reads_to_remove:
                putative_cluster_b.remove(read)
        return putative_cluster_a, putative_cluster_b

    @property
    def orientation_vector(self):
        """Return orientation of all reads with 'BD' tag."""
        return [('R', i) if r.is_reverse else ('F', i) for i, r in enumerate(self) if not r.has_tag('AD')]

    @property
    def read_index(self):
        """Index of read names in cluster."""
        return set([r.query_name for r in self if not r.has_tag('AC')])  # AC is assembled contig

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
            self._clustertag = TagCluster(self, shm_dir=self.shm_dir)
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
        other_clustertag = TagCluster(other_cluster, shm_dir=self.shm_dir)
        # TODO: ... the should be no additional polarity switch if clusters are to be joined?
        # Second check ... are three_p and five_p of cluster overlapping?
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
        # TODO: Refactor this to a common function for left-left and right-right assembly
        if self.clustertag.right_sequences and other_clustertag.right_sequences:
            # A cluster that provides support for a 5p insertion will have the reads always annotated as left sequences.
            # That's a bit confusing, since the mates are to the right of the cluster ... but that's how it is.
            if (other_clustertag.three_p_breakpoint - self.clustertag.three_p_breakpoint) < max_distance:
                # We don't want clusters to be spaced too far away. Not sure if this is really a problem in practice.
                if Cap3Assembly.sequences_contribute_to_same_contig(self.clustertag.right_sequences, other_clustertag.right_sequences):
                    self._can_join_d = {self.hash: other_cluster.hash}
                    return True
        self_switches = self.orientation_switches
        other_switches = other_cluster.orientation_switches
        if len(self_switches) == 1 and len(other_switches) and self_switches != other_switches:
            # Merge a cluster that is split by an insertion and not connected via split reads
            if self_switches[0][0] == 'F' and other_switches[0][0] == 'R':
                self._can_join_d = {self.hash: other_cluster.hash}
                return True
        # We know this cluster (self) cannot be joined with other_cluster, so we cache this result,
        # Since we may ask this question multiple times when joining the clusters.
        self._cannot_join_d = {self.hash: other_cluster.hash}
        return False

    @property
    def left_support(self):
        """Number of supporting reads on the left side of cluster."""
        return self.clustertag.left_sequence_count

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
        # We just touch the contigs in threading mode,
        # so the actual assembly is being triggered.
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

    @property
    def start_corrected(self):
        """
        Extend start site.

        May be extended by maximum insert size length minus the most distant 5' mate if no exact TSD exists.
        """
        if not self.clustertag.tsd.three_p_support and len(self.orientation_switches) < 2:
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
            three_p_reads = self.left_mate_support.values()
            if three_p_reads:
                start_position = min([r.reference_start for r in three_p_reads])
                end_position = max([r.reference_end for r in three_p_reads])
                if end_position - start_position < self.max_proper_size:
                    return start_position + self.max_proper_size
        return self.clustertag.three_p_breakpoint

    @cached_property
    def _start_and_end(self):
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
    min_start = start - 500
    if min_start <= 0:
        min_start = 1
    end = last_chunk[2]
    max_end = end + 500
    try:
        f = pysam.AlignmentFile(input_path)
        try:
            reads = f.fetch(chromosome, min_start, max_end)
        except Exception:
            pysam.index(input_path)
            f = pysam.AlignmentFile(input_path)
            reads = f.fetch(chromosome, min_start, max_end)
        for r in reads:
            if not r.is_duplicate and r.is_proper_pair and r.mapq > 0:
                add_to_clusters(chunk, r, result)
        f.close()
    except ValueError:
        logging.warn("Encountered ValueError on chromosome %s for start %s and end %s of chunks %s" % (chromosome, start, end, chunk))
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
