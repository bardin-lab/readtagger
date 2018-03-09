from collections import defaultdict


AMBIGUOUS_IUPAC = {'A': 'A',
                   'C': 'C',
                   'G': 'G',
                   'T': 'T',
                   'N': 'N',
                   'AG': 'R',
                   'CT': 'Y',
                   'CG': 'S',
                   'AT': 'W',
                   'GT': 'K',
                   'AC': 'M',
                   'CGT': 'B',
                   'AGT': 'D',
                   'ACT': 'H',
                   'ACG': 'V',
                   'ACGT': 'N'}

IUAPC_AMBIGUOUS = {v: k for k, v in AMBIGUOUS_IUPAC.items()}


def dumb_consensus(string_list, left_align=True):
    """
    Get the dumb consensus from a list of nucelotide strings.

    :param string_list:  list of sequences for which to get the consensus
    :param left_align:   if True get consensus from left to right, otherwise from right to left
    :return:

    >>> dumb_consensus(['ATGC'])
    'ATGC'
    >>> dumb_consensus(['ATGC', 'ATGC', 'ATGG'])
    'ATGC'
    >>> dumb_consensus(['ATGC', 'ATGC', 'AGTT'])
    'ATGC'
    >>> dumb_consensus(['ATGC', 'ATGT'])
    'ATGY'
    >>> dumb_consensus(['ATNC', 'ATNT'])
    'ATNY'
    >>> dumb_consensus(['ATGC', 'ATGA', 'ATGT', 'ATGG'])
    'ATGN'
    >>> dumb_consensus(['ATGC', 'ATGA', 'ATGT', 'ATGN'])
    'ATGH'
    >>> dumb_consensus(['ATGCA', 'ATGC'])
    'ATGCA'
    >>> dumb_consensus(['AAAATGC', 'ATTGC', 'TATGC', 'AATGC'], left_align=False)
    'AAAATGC'
    """
    new_consensus = []
    assert isinstance(string_list, list)
    max_len = max(len(l) for l in string_list)
    if not left_align:
        string_list = right_align_strings(string_list=string_list)
    for p in range(max_len):
        p_nt_occurrence = defaultdict(int)
        for l in string_list:
            try:
                p_nt_occurrence[l[p]] += 1
            except IndexError:  # we are in one of the shorter strings
                continue
        most_common = max(p_nt_occurrence, key=lambda key: p_nt_occurrence[key])
        ties = {key for key, value in p_nt_occurrence.items() if value == p_nt_occurrence[most_common]}
        if len(ties) > 1 and 'N' in ties:
            ties -= {'N'}
        # Lookup IUPAC code for ambiguous bases
        new_consensus.append(AMBIGUOUS_IUPAC["".join(sorted(ties))])
    if left_align:
        return "".join(new_consensus)
    else:
        return "".join(new_consensus[::-1])


def right_align_strings(string_list):
    """Invert the sequence of a nested list."""
    right_aligned_string_list = []
    for l in string_list:
        right_aligned_string_list.append(l[::-1])
    return right_aligned_string_list
