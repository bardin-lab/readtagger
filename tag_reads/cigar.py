import re

pattern = re.compile('([MIDNSHPX=])')
#  Op BAM Description
#  ------------------
#  M 0 alignment match (can be a sequence match or mismatch)
#  I 1 insertion to the reference
#  D 2 deletion from the reference
#  N 3 skipped region from the reference
#  S 4 soft clipping (clipped sequences present in SEQ)
#  H 5 hard clipping (clipped sequences NOT present in SEQ)
#  P 6 padding (silent deletion from padded reference)
#  = 7 sequence match
#  X 8 sequence mismatch

# M includes both sequence match and mismatch (that's why = and X exist)
# Seems BWA is not setting those (same for N, which is being set by tophat/hisat/star for introns)


def cigar_to_tuple(cigar, reverse=False):
    """
    Turns a cigar into a list of tuples by splitting at any of the valid CIGAR letters
    (MIDNSHPX)
    :param cigar: string
    :return: list of tuples, with first integer spcifying the length of the operation,
              and the string specifying the operation
    >>> cigar_to_tuple('3M1I44M1D7M1I32M37S')
    [(3, 'M'), (1, 'I'), (44, 'M'), (1, 'D'), (7, 'M'), (1, 'I'), (32, 'M'), (37, 'S')]
    >>> cigar_to_tuple('3M1I44M1D7M1I32M37S', reverse=True)
    [(37, 'S'), (32, 'M'), (1, 'I'), (7, 'M'), (1, 'D'), (44, 'M'), (1, 'I'), (3, 'M')]
    """
    result = pattern.split(cigar)[:-1]
    result = zip(map(int, result[0::2]), result[1::2])
    if reverse:
        result.reverse()
    return result


def cigar_tuple_to_cigar_length(cigar):
    """
    Replace consecutive tags with coordinates in the read
    >>> cigar = [(3, 'M'), (1, 'I'), (44, 'M'), (1, 'D'), (7, 'M'), (1, 'I'), (32, 'M'), (37, 'S')]
    >>> cigar_tuple_to_cigar_length(cigar)
    [((0, 3), 'M'), ((3, 4), 'I'), ((4, 48), 'M'), ((48, 49), 'D'), ((49, 56), 'M'), ((56, 57), 'I'), ((57, 89), 'M'), ((89, 126), 'S')]
    """
    cigar_length_tuples = []
    start = 0
    for (i, m) in cigar:
        end = start + i
        cigar_length_tuples.append(((start, end), m))
        start = end
    return cigar_length_tuples


def get_cigar_lengths(cigar, reverse=False):
    return cigar_tuple_to_cigar_length(cigar_to_tuple(cigar, reverse=reverse))


def alternative_alignment_cigar_is_better(current_cigar, alternative_cigar, same_orientation):
    """
    Try to judge whether a read should be marked as having an alternative alignment
    tag. The read with the characteristics below should not be marked as aligning
    to a transposon, because:
      - This read has 10 nt softclipped on the left side, and the alternative transposon-aligned
        read is antisense with a cigar of 43M1D7M1I32M42S. Because the read is mapping antisense to the TE,
        we need to revert the order of the CIGAR string if we want to compare the CIGAR.
        The TE alignment therefor starts with 42 softclipped bases and should be discarded, because
        the suboptimal alignment on the reference genome starts with 10 soft-clipped bases. The alternative alignment
        does not provide a better explanation for the clipping we see.
      - The mate is also marked as aligning to a TE (in antisense orientation) and starts with 6 softclipped
        nucleotides, followed by 39 matches and 80 softclipped nucleotides. This is the same 6 nucleotides of clipping
        between this alignment and the alternative alignment, but the current aligned segment explains the remaining
        matches nucleotides much better, so the alternative tag should also be removed.

     In [9]: r.query_name
     Out[9]: 'HWI-D00405:129:C6KNAANXX:4:1110:14319:63440'
     In [7]: r.cigarstring
     Out[7]: '10S115M'
     In [8]: r.is_reverse
     Out[8]: False
     In [9]: r.tags
     Out[9]:
     [('AD', 'R:FBti0060251_Dmel\\1360_P,POS:666,CIGAR:43M1D7M1I32M42S,S:AS'),
     ('AR', 'FBti0060251_Dmel\\1360_P'),
     ('BD', 'R:FBti0060251_Dmel\\1360_P,POS:665,CIGAR:80S39M6S,S:S'),
     ('BR', 'FBti0060251_Dmel\\1360_P'),
     ('MC', '6S119M'),
     :param cigar: cigar string of current read
     :param alternative_cigar: alternative cigar
     :return: bool

     >>> ad_args = dict(current_cigar='10S115M', alternative_cigar='43M1D7M1I32M42S', same_orientation=False)
     >>> alternative_alignment_cigar_is_better(**ad_args)
     False
     >>> bd_args = dict(current_cigar='6S119M', alternative_cigar='80S39M6S', same_orientation=False)
     >>> alternative_alignment_cigar_is_better(**bd_args)
     False
     >>> true_alternatives = [('94M31S', '32M93S', False), ('95M30S', '90S35M', True), ('95M30S', '90S35M', True), ('92M33S', '87S38M', True), ('92M33S', '87S38M', True), ('65M60S', '66M59S', False), ('31S94M', '32M93S', True)]
     >>> [alternative_alignment_cigar_is_better(*a) for a in true_alternatives]
     [True, True, True, True, True, True, True]
    """
    # We start by discarding alignments that are not clipped
    if ('S' or 'H') not in current_cigar:
        return False
    if ('S' or 'H') not in alternative_cigar:
        # I don't think there is a case in which we would want to discard an alternative cigar that is better than the
        # current cigar (no clipping in alternative read, while current cigar is clipped) UNLESS the alternative cigar
        # has large deletions or insertions ... maybe we should at least check that the amount of matched bases is higher
        #  TODO: verify this.
        return True
    cigar_lengths_current_read = get_cigar_lengths(current_cigar)
    split_regions_current_read = [set(range(*reg)) for (reg, operation) in cigar_lengths_current_read if operation in {'S', 'H'}]  # These are the split regions wrt the current BAM file
    cigar_lengths_alternative = get_cigar_lengths(alternative_cigar, reverse=not same_orientation)
    matched_regions_alternative = [set(range(*reg)) for (reg, operation) in cigar_lengths_alternative if operation == 'M']
    for reg in split_regions_current_read:
        for altreg in matched_regions_alternative:
            if reg & altreg:
                # clipped region in current read overlaps with aligned region in alternative alignment.
                # this does includes situations like cigar_reg = 3S117M and cigar_altreg = 19M101S
                # Not sure if we want to further filter those out ...  .
                # Ideally the cigar_altreg M operation would cover the length of the softclipped portion
                # (or slightly more.) That's frequently what I see for True insertions (chr3R:13,762,417-13,767,183, HUM4)
                # The only case that's not true is if there are minor insertions or deletions in the transposon aligned
                # read. So perhaps I should stitch together Matched regions that are separated by I or D events
                # For I events the pseudo M should get longer by the length of the I event, and shorter by the length of the D
                # event. Will check this later.
                if len(altreg) >= len(reg):
                    # The alternative match overlaps the clipped region and is equal to or longer than the clipped region
                    return True
    return False






