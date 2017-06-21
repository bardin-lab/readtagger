from collections import Counter
from cached_property import cached_property
from .tags import Tag

MAX_TSD_SIZE = 50
HARD_CLIP = 5


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

    @property
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
        return cigar[0][0] == HARD_CLIP

    @staticmethod
    def _hard_clip_right(read):
        cigar = Tag.from_read(read).cigar
        return cigar[-1][0] == HARD_CLIP

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
