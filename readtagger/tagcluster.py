"""Organize a cluster of interest."""
from collections import Counter
from .cigar import cigar_tuple_to_cigar_length
from .cigar import stitch_matched_regions
from .tags import Tag


class TagCluster(object):
    """Take a cluster of reads and figure out how to organize clustered reads."""

    def __init__(self, cluster):
        """Cluster is an iterable of pysam.AlignedSegement objects."""
        self.cluster = cluster
        self.sequences_of_interest = self.extract_seqs_of_interest()

    def find_break(self):
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
        :return:
        """
        # Start by finding reads with informative splits
        pass

    def extract_seqs_of_interest(self):
        """Extract sequences for further processing, like assembly or MSA."""
        # Mate sequences for reads whose mate has an alternative alignment
        sequences = {"%s_m" % (r.query_name): r.get_tag('MS') for r in self.cluster if r.has_tag('BD')}
        splits = [r for r in self.cluster if r.has_tag('AD')]
        split_tags = [Tag.from_tag_str(r.get_tag('AD')) for r in splits]
        cigar_tuples = [tag.cigar[::-1] if r.is_reverse != tag.is_reverse else tag.cigar for tag, r in zip(split_tags, splits)]
        cigar_lengths = [stitch_matched_regions(cigar_tuple_to_cigar_length(c)) for c in cigar_tuples]
        split_sequences = {r.query_name: r.query_sequence[start:end] for c, r in zip(cigar_lengths, splits) for ((start, end), op) in c if op == 0}
        sequences.update(split_sequences)
        return sequences


class TargetSiteDuplication(object):
    """Collect Target Site Duplication methods and data."""

    def __init__(self, cluster):
        self.cluster = cluster
        self.split_ads = [r for r in self.cluster if r.has_tag('AD')]
        self.three_p = self.find_three_p()
        self.five_p = self.find_five_p()

    def find_five_p(self):
        """
        Find the five prime genomic position of a TSD.

        It is frequently further 3p than the three prime of the same TSD.
        See the bewlow situation, read marked with * supports a TSD five prime, X indicates clipping.

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
        reference_ends = sorted([r.reference_end for r in self.split_ads])
        min_five_p_aligned_end = min(reference_ends)
        most_common_ends = Counter(reference_ends).most_common()
        most_common_end_pos_occurence = most_common_ends[0][1]
        best_ends = [t[0] for t in most_common_ends if t[1] == most_common_end_pos_occurence]
        if min_five_p_aligned_end not in best_ends:
            # We could have a rare misaligned softclipped read (that happens if there is a mismatch close to the acutal clipped region),
            # so we check if the next position (in the genomic 3p direction) may be better (= the most frequent site).
            for p in reference_ends:
                if p == min_five_p_aligned_end:
                    continue
                if p in best_ends:
                    min_five_p_aligned_end = p
                    break
                else:
                    break
        return min_five_p_aligned_end

    def find_three_p(self):
        """
        Find the five prime genomic position of a TSD.

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
        starting_positions = sorted([r.pos for r in self.split_ads])
        max_starting_position = max(starting_positions)
        most_common_starts = Counter(starting_positions).most_common()
        most_common_occurence = most_common_starts[0][1]
        best_starts = [t[0] for t in most_common_starts if t[1] == most_common_occurence]
        if max_starting_position not in best_starts:
            # Normally the 3p soft clipped position for a TE insertion should be
            # - overrepresented among the starting positions of clipped reads
            # - the rightmost starting position
            # If that's not the case for certain reads (e.g because there is a mismatch closeby), check
            # if perhaps the next position is the mode and use this.
            for p in starting_positions[::-1]:
                if p == max_starting_position:
                    continue
                if p in best_starts:
                    max_starting_position = p
                    break
                else:
                    break
        return max_starting_position

    @property
    def three_p_support(self):
        if not hasattr(self, '_three_p_support'):
            self._three_p_support = [r.query_name for r in self.split_ads if r.pos == self.three_p]
        return self._three_p_support

    @property
    def five_p_support(self):
        if not hasattr(self, '_five_p_support'):
            self._five_p_support = [r.query_name for r in self.split_ads if r.reference_end == self.five_p]
        return self._five_p_support

    @property
    def unassigned_support(self):
        if not hasattr(self, '_unassigned_support'):
            self._unassigned_support = [r.query_name for r in self.split_ads if r.query_name not in self.five_p_support + self.three_p_support]
        return self._unassigned_support
