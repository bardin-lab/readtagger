from itertools import groupby

CODE2CIGAR = "MIDNSHP=XB"
CIGAR2CODE = dict([y, x] for x, y in enumerate(CODE2CIGAR))


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


def cigartuples_to_cigarstring(cigartuple):
    """
    Convert cigartuple to cigarstring.

    :param cigartuple:
    :rtype: str
    >>> cigartuples_to_cigarstring([(0, 91), (4, 34)])
    '91M34S'
    """
    return "".join(["%s%s" % (l, CODE2CIGAR[op]) for op, l in cigartuple])


def cigar_split(cigarstring):
    """
    Group cigarstring.

    CIGAR grouping function modified from https://github.com/brentp/bwa-meth
    and https://github.com/mdshw5/simplesam/blob/master/simplesam.py .

    :param cigarstring: Cigarstring to iterate over
    :rtype generator
    """
    cig_iter = groupby(cigarstring, lambda c: c.isdigit())
    for _, n in cig_iter:
        yield ("".join(n)), "".join(next(cig_iter)[1], )


def cigar_to_tuple(cigar):
    """
    Turn a cigar into a list of tuples.

    :param cigar: string
    :return: list of tuples, with first integer spcifying the length of the operation,
              and the string specifying the operation
    >>> cigar_to_tuple('3M1I44M1D7M1I32M37S')
    [(0, 3), (1, 1), (0, 44), (2, 1), (0, 7), (1, 1), (0, 32), (4, 37)]
    >>> c = cigar_to_tuple('3M1I44M1D7M1I32M37S')
    >>> c.reverse()
    >>> c
    [(4, 37), (0, 32), (1, 1), (0, 7), (2, 1), (0, 44), (1, 1), (0, 3)]
    """
    return [(CIGAR2CODE[op], int(l)) for (l, op) in cigar_split(cigar)]


def cigar_tuple_to_cigar_length(cigar):
    """
    Replace consecutive tags with coordinates in the read.

    >>> cigar = [(0, 3), (1, 1), (0, 44), (3, 1), (0, 7), (1, 1), (0, 32), (4, 37)]
    >>> cigar_tuple_to_cigar_length(cigar)
    [((0, 3), 0), ((3, 4), 1), ((4, 48), 0), ((48, 49), 3), ((49, 56), 0), ((56, 57), 1), ((57, 89), 0), ((89, 126), 4)]
    """
    cigar_length_tuples = []
    start = 0
    for (m, i) in cigar:
        end = start + i
        cigar_length_tuples.append(((start, end), m))
        start = end
    return cigar_length_tuples


def alternative_alignment_cigar_is_better(current_cigar, alternative_cigar, same_orientation):
    """
    Judge whether a read should be marked as having an alternative alignment tag.

    The read with the characteristics below should not be marked as aligning
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
     [('AD', 'R:FBti0060251_Dmel1360_P,POS:666,CIGAR:43M1D7M1I32M42S,S:AS'),
     ('AR', 'FBti0060251_Dmel1360_P'),
     ('BD', 'R:FBti0060251_Dmel1360_P,POS:665,CIGAR:80S39M6S,S:S'),
     ('BR', 'FBti0060251_Dmel1360_P'),
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
     >>> cd_args = dict(current_cigar='6S119M', alternative_cigar='6S119M', same_orientation=True)
     >>> alternative_alignment_cigar_is_better(**cd_args)
     False
          >>> cd_args = dict(current_cigar='6S119M', alternative_cigar='119M6S', same_orientation=False)
     >>> alternative_alignment_cigar_is_better(**cd_args)
     False
     >>> true_alternatives = [('94M31S', '32M93S', False), ('95M30S', '90S35M', True), ('95M30S', '90S35M', True)]
     >>> [alternative_alignment_cigar_is_better(*a) for a in true_alternatives]
     [True, True, True]
    """
    # We start by discarding alignments that are not clipped
    if isinstance(current_cigar, str):
        current_cigar = cigar_to_tuple(current_cigar)
    if isinstance(alternative_cigar, str):
        alternative_cigar = cigar_to_tuple(alternative_cigar)
    if not same_orientation:
        alternative_cigar = list(alternative_cigar)
        alternative_cigar.reverse()
    if not any([op for (op, l) in current_cigar if op in {4, 5}]):
        # No soft or hard-clipping in current cigar. This is a perfect aignment, alternative cigar can't be better
        return False
    if not any([op for (op, l) in alternative_cigar if op in {4, 5}]):
        # No clipping in alternative cigar, while current cigar is clipped.
        # I don't think there is a case in which we would want to discard an alternative cigar that is better than the
        # current cigar UNLESS the alternative cigar has large deletions or insertions.
        # Maybe we should at least check that the amount of matched bases is higher
        #  TODO: verify this.
        return True
    cigar_lengths_current_read = cigar_tuple_to_cigar_length(current_cigar)
    # Next we get a range representing the soft and hardlclipped regions in the current read
    split_regions_current_read = [set(range(*reg)) for (reg, operation) in cigar_lengths_current_read if operation in {4, 5}]
    cigar_lengths_alternative = cigar_tuple_to_cigar_length(alternative_cigar)
    # Next are the matched regions in the alternative cigar.
    matched_regions_alternative = [set(range(*reg)) for (reg, operation) in cigar_lengths_alternative if operation == 0]
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
