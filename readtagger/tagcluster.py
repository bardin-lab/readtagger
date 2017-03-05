"""Organize a cluster of interest."""
from collections import Counter
import warnings
from cached_property import cached_property
from .cap3 import Cap3Assembly
from .tags import Tag

MAX_TSD_SIZE = 50


class TagCluster(object):
    """
    Take a cluster of reads and figure out how to organize clustered reads.

    First attempts to find the Target Site Duplication (TSD). Then tries to get the breakpoint,
    and if we get the breakpoint we can extract the sequences left and right of the insertion and
    all sequences pointing into the insertion to assemble the putative inserted sequence, for the left and right side.
    """

    def __init__(self, cluster):
        """Cluster is an iterable of pysam.AlignedSegment objects."""
        self.cluster = cluster
        self.tsd = TargetSiteDuplication(self.cluster)
        self.five_p_breakpoint, self.three_p_breakpoint = self.find_breakpoint()
        # TODO: Implement looking up what TE the insert(s) best match to
        # TODO: Sum read support:
        #   - how many reads support TSD?
        #   - how many fragments cover any part of the insertion?

    @cached_property
    def left_insert(self):
        """Return insert sequence as assembled from the left side."""
        if self.left_sequences:
            if not hasattr(self, '_left_seq_cap3'):
                self._left_seq_cap3 = Cap3Assembly(self.left_sequences)
            return self._left_seq_cap3
        else:
            return None

    @cached_property
    def right_insert(self):
        """Return insert sequence as assembled from the right side."""
        if self.right_sequences:
            if not hasattr(self, '_right_seq_cap3'):
                self._right_seq_cap3 = Cap3Assembly(self.right_sequences)
            return self._right_seq_cap3
        else:
            return None

    @cached_property
    def joint_insert(self):
        """Return joint insert sequence."""
        if self.right_sequences and self.left_sequences:
            return Cap3Assembly.join_assemblies([self.left_insert, self.right_insert])
        elif self.right_sequences:
            return self.right_insert
        elif self.left_sequences:
            return self.left_insert
        else:
            return None

    def find_breakpoint(self):
        """
        Find the breakpoint of a potential insertion.

        We know that informative reads on the left of a insertion that don't overlap the insertion (i.e that are not split)
        must always be sense. Inversely, reads that are on the right of an insertion must always be antisense.
        Reads that are split should have their split region pointing into the insertion.

        Genome:                     |012345678901234567890
        Insertion:                  |----------vv----------
        Read left, Mate in TE:      |>>>>>>>>>>
        Read left, split in TE:     |    >>>>>>>>XX
        Target Site Duplication TSD:|          TT
        Read right, split in TE:    |        XX<<<<<<<<
        Read right, Mate in TE:     |           <<<<<<<<<<

        We have multiple options to find the insertion breakpoint.
         - We can identify the Target Site Duplication.
         - We can determine the closest read pairs before the sense of the aligned mates switches.
         - We can cluster the split and informative mates and then infer from the 2 (or more ...) clusters that are generated
           where the insertion should be placed.
        If we rely on TSD detection we may loose some information we can get from mate pairs.
        The last option will work for long read sequencing while potentially keeping the paired end information.
        """
        # Start by finding reads with informative splits
        if not self.tsd.is_valid:
            # Need to determine the break more closely.
            # If 3' or 5' of breakpoint has been identified -- use that.
            # TODO: Could eventually be improved by looking for softclipped positions without AD,
            # but then I would need to look at reads that are not marked to be in the cluster.
            # I could look at the BD tagged reads though and scan for softclipping ... .
            five_p, three_p = self.tsd.five_p, self.tsd.three_p
            if three_p is None:
                three_p = self.infer_three_p_from_mates()
            if five_p is None:
                five_p = self.infer_five_p_from_mates()
            if self.tsd.five_p and self.tsd.three_p and self.tsd.three_p < self.tsd.five_p:
                # An invalid TSD was found, probably because the inferred TSD is too long.
                # This could be two close-by insertions that each have support from one side,
                # and/or --less likely-- a pre-existing duplication of that sequence?
                # Is that even a single insertion? Perhaps check that evidence
                # points to the same TE? Should I split the cluster? Maybe I should just dump
                # the cluster reads as a BAM file and inspect them to see what should be done.
                warn = "Found a cluster with 5p and 3p evidence for TSD, but reads are spaced too far apart.\n"
                warn += "The cluster coordinates are tid: %s, start:%s, end: %s" % (self.cluster[0].tid, self.cluster[0].pos, self.cluster[-1].reference_end)
                warnings.warn(warn)
            return five_p, three_p
        return self.tsd.five_p, self.tsd.three_p

    def infer_five_p_from_mates(self):
        """Return rightmost reference end for sense reads with BD tag."""
        five_p_reads = [r.reference_end for r in self.cluster if not r.is_reverse and r.has_tag('BD')]
        if five_p_reads:
            return max(five_p_reads)
        else:
            return None

    def infer_three_p_from_mates(self):
        """Return leftmost reference start for antisense reads with BD tag."""
        three_p_reads = [r.pos for r in self.cluster if r.is_reverse and r.has_tag('BD')]
        if three_p_reads:
            return min(three_p_reads)
        else:
            return None

    @cached_property
    def left_sequences(self):
        """
        Find reads left of a breakpoint.

        These reads need to be sense oriented if they have a BD tag, and reads with an
        AD tag should support that particular TSD end (5p for left reads).
        """
        # TODO: extend this to also return quality values
        left_sequences = {}
        for r in self.cluster:
            if r.has_tag('BD'):
                if not r.is_reverse:
                    if r.is_read1:
                        qname = "%s.1" % r.query_name
                    else:
                        qname = "%s.2" % r.query_name
                    left_sequences[qname] = r.get_tag('MS')
            if r.has_tag('AD') and r.query_name in self.tsd.five_p_support:
                if r.reference_end == self.tsd.five_p or r.pos == self.tsd.five_p:
                    left_sequences[r.query_name] = r.query_sequence
        return left_sequences

    @cached_property
    def right_sequences(self):
        """
        Find reads right of a breakpoint.

        These reads need to be antisense oriented if they have a BD tag, and reads with an
        AD tag should support that particular TSD end (3p for right reads).
        """
        # TODO: extend this to also return quality values
        right_sequences = {}
        for r in self.cluster:
            if r.has_tag('BD'):
                # Exclude reads that have both AD and BD when determining reads that support the right side,
                # because if they overlap they will already count by their AD Tag below.
                if r.is_reverse:
                    if r.is_read1:
                        qname = "%s.1" % r.query_name
                    else:
                        qname = "%s.2" % r.query_name
                    right_sequences[qname] = r.get_tag('MS')
            if r.has_tag('AD') and r.query_name in self.tsd.three_p_support:
                if r.reference_end == self.tsd.three_p or r.pos == self.tsd.three_p:
                    right_sequences[r.query_name] = r.query_sequence
        return right_sequences


class TargetSiteDuplication(object):
    """Collect Target Site Duplication methods and data."""

    def __init__(self, cluster, include_duplicates=False):
        """
        Return TargetSiteDuplication object for reads in cluster.

        >>> from test.helpers import MockAlignedSegment as R
        >>> # 2 nt 3' clipping
        >>> r1 = R(query_name='r1', reference_end=10, pos=2, query_alignment_start=0, query_alignment_end=10, query_length=12)
        >>> r2 = R(query_name='r2', reference_end=12, pos=2, query_alignment_start=0, query_alignment_end=12, query_length=15)
        >>> r3 = R(query_name='r3', reference_end=12, pos=2, query_alignment_start=0, query_alignment_end=12, query_length=15)
        >>> # 3 nt 3' clipping
        >>> r4 = R(query_name='r4', reference_end=12, pos=2, query_alignment_start=0, query_alignment_end=12, query_length=15)
        >>> tsd = TargetSiteDuplication(cluster = [r1, r2, r3, r4], include_duplicates=True)
        >>> assert tsd.three_p is None
        >>> tsd.five_p
        12
        >>> tsd.is_valid
        False
        >>> # 10 nt 5' clipping
        >>> r5 = R(query_name='r5', reference_end=20, pos=10, query_alignment_start=10, query_alignment_end=20, query_length=20)
        >>> r6 = R(query_name='r6', reference_end=20, pos=10, query_alignment_start=10, query_alignment_end=20, query_length=20)
        >>> r7 = R(query_name='r7', reference_end=20, pos=10, query_alignment_start=10, query_alignment_end=20, query_length=20)
        >>> # 12 nt 5' clipping
        >>> r8 = R(query_name='r8', reference_end=20, pos=12, query_alignment_start=12, query_alignment_end=20, query_length=20)
        >>> tsd = TargetSiteDuplication(cluster = [r5, r6, r7, r8], include_duplicates=True)
        >>> assert tsd.five_p is None
        >>> tsd.three_p
        10
        >>> tsd.is_valid
        False
        >>> tsd = TargetSiteDuplication(cluster = [r1, r2, r3, r4, r5, r6, r7, r8], include_duplicates=True)
        >>> tsd.is_valid
        True
        >>> tsd.five_p_clip_length
        3
        >>> tsd.three_p_clip_length
        10
        >>> tsd.three_p_support
        ['r5', 'r6', 'r7']
        >>> tsd.five_p_support
        ['r2', 'r3', 'r4']
        >>> tsd.unassigned_support
        ['r1', 'r8']
        """
        self.cluster = cluster
        if include_duplicates:
            self.split_ads = [r for r in self.cluster if r.has_tag('AD')]
        else:
            self.split_ads = [r for r in self.cluster if r.has_tag('AD') and not r.is_duplicate]
        self.three_p = self.find_three_p()
        self.five_p = self.find_five_p()

    @cached_property
    def is_valid(self):
        """
        Return True if Target Site Duplication is valid.

        A TSD is valid if:
          - A three prime and five prime site have been found in an instance of this class.
            This implies that five prime and three prime are defined by split reads
          - The three prime must be left to the five prime
          - The length of the TSD must be less then 51 nucleotides
        """
        # Super arbitrary, but I guess this is necessary
        if self.three_p and self.five_p and self.three_p != self.five_p and self.three_p < self.five_p and (abs(self.five_p - self.three_p) <= MAX_TSD_SIZE):
            return True
        else:
            return False

    def find_five_p(self):
        """
        Find the five prime genomic position of a TSD.

        It is frequently further 3p than the three prime of the same TSD.
        See the below situation, read marked with * supports a TSD five prime, X indicates clipping.

          Genome:                     |012345678901234567890
          Insertion:                  |----------vv----------
          Read left, Mate in TE:      |>>>>>>>>>>
        * Read left, split in TE:     |    >>>>>>>>XX
          Target Site Duplication TSD:|          TT
          Read right, split in TE:    |        XX<<<<<<<<
          Read right, Mate in TE:     |           <<<<<<<<<<

        The five prime of a TSD is characterized by having
          - the rightmost alignment end (of reads that support an insertion with the AD tag)
          - the most frequent alignment end (of reads that support an insertion with the AD tag)
        """
        if not self.sorted_split_end_positions:
            return None
        min_five_p_aligned_end = min(self.sorted_split_end_positions)
        most_common_ends = Counter(self.sorted_split_end_positions).most_common()
        most_common_end_pos_occurence = most_common_ends[0][1]
        best_ends = [t[0] for t in most_common_ends if t[1] == most_common_end_pos_occurence]
        if min_five_p_aligned_end not in best_ends:
            # We could have a rare misaligned softclipped read (that happens if there is a mismatch close to the acutal clipped region),
            # so we check if the next position (in the genomic 3p direction) may be better (= the most frequent site).
            for p in self.sorted_split_end_positions:
                if p == min_five_p_aligned_end:
                    continue
                if p in best_ends:
                    return p
                else:
                    break
        return min_five_p_aligned_end

    def find_three_p(self):
        """
        Find the three prime genomic position of a TSD.

        It is frequently further 5p than the five prime of the same TSD.
        See the below situation, read marked with * supports a TSD five prime, X indicates clipping.

          Genome:                     |012345678901234567890
          Insertion:                  |----------vv----------
          Read left, Mate in TE:      |>>>>>>>>>>
          Read left, split in TE:     |    >>>>>>>>XX
          Target Site Duplication TSD:|          TT
        * Read right, split in TE:    |        XX<<<<<<<<
          Read right, Mate in TE:     |           <<<<<<<<<<

        The three prime of a TSD is characterized by having
          - the leftmost alignment start (of reads that support an insertion with the AD tag)
          - the most frequent alignment start (of reads that support an insertion with the AD tag)
        """
        if not self.sorted_split_start_positions:
            return None
        max_starting_position = max(self.sorted_split_start_positions)
        most_common_starts = Counter(self.sorted_split_start_positions).most_common()
        # Counter().most_common() returns a list of tuples,
        # where tuples are ordered from most to least common,
        # and where the first item in a tuple is the value
        # and the second is the occurence.
        most_common_occurence = most_common_starts[0][1]
        best_starts = [t[0] for t in most_common_starts if t[1] == most_common_occurence]
        if max_starting_position not in best_starts:
            # Normally the 3p soft clipped position for a TE insertion should be
            # - overrepresented among the starting positions of clipped reads
            # - the rightmost starting position
            # If that's not the case for certain reads (e.g the clip has been wrongly extended because there is an adjacent mismatch),
            # check if perhaps the next position is the mode and use this.
            for p in self.sorted_split_start_positions[::-1]:
                if p == max_starting_position:
                    continue
                if p in best_starts:
                    return p
                else:
                    break
        return max_starting_position

    @staticmethod
    def _hard_clip_left(read):
        cigar = Tag.from_read(read).cigar
        if cigar[0][0] == 5:
            return True
        else:
            return False

    @staticmethod
    def _hard_clip_right(read):
        cigar = Tag.from_read(read).cigar
        if cigar[-1][0] == 5:
            return True
        else:
            return False

    @cached_property
    def sorted_split_start_positions(self):
        """Return sorted start positions of split reads."""
        return sorted([r.pos for r in self.split_ads if r.query_alignment_start != 0 or self._hard_clip_left(r)])  # To only get relevant split start positions

    @cached_property
    def sorted_split_end_positions(self):
        """Return sorted end positions of split reads."""
        return sorted([r.reference_end for r in self.split_ads if r.query_alignment_end != r.query_length or self._hard_clip_right(r)])

    @cached_property
    def three_p_clip_length(self):
        """
        Return the longest clipped region that extends the three_p breakpoint.

        5p Breakpoint:           v
        3p Breakpoint:               v
        Genome:         |012345678901234567890
        5p split read:  |  >>>>>>XXXX
        3p split read:  |          XXX>>>>>>>

        In this case return 3 for the three_p breakpoint.
        """
        # TODO: This doesn't work for hardclipped reads. Should use left or right clip
        clip_length_at_three_p = [r.query_alignment_start for r in self.split_ads if r.pos == self.three_p]
        if not clip_length_at_three_p:
            return 0
        return max(clip_length_at_three_p)

    @cached_property
    def five_p_clip_length(self):
        """
        Return the longest clipped region that extends the five_p breakpoint.

        5p Breakpoint:           v
        3p Breakpoint:               v
        Genome:         |012345678901234567890
        5p split read:  |  >>>>>>XXXX
        3p split read:  |          XXX>>>>>>>

        In this case return 4 for the 5p breakpoint.
        """
        # TODO: This doesn't work for hardclipped reads. Should use left or right clip
        clip_length_at_five_p = [r.query_length - r.query_alignment_end for r in self.split_ads if r.reference_end == self.five_p]
        if not clip_length_at_five_p:
            return 0
        return max(clip_length_at_five_p)

    @cached_property
    def three_p_support(self):
        """Return list of Reads that support the inferred three prime position for this TSD."""
        if not hasattr(self, '_three_p_support'):
            self._three_p_support = [r.query_name for r in self.split_ads if r.pos == self.three_p]
        return self._three_p_support

    @cached_property
    def five_p_support(self):
        """Return list of Reads that support the inferred five prime position for this TSD."""
        return [r.query_name for r in self.split_ads if r.reference_end == self.five_p]

    @cached_property
    def unassigned_support(self):
        """Return list of Reads that were not starting or ending at the three prime or five prime of this TSD."""
        return [r.query_name for r in self.split_ads if r.query_name not in self.five_p_support + self.three_p_support]
