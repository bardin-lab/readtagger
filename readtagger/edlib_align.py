import re
try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

from edlib import align


COMPLEMENTARY_SEQUENCES = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}


def multiple_sequences_overlap(queries, targets, min_match=30):
    """Check if any query in queries overlaps with any target in targets."""
    for seq in queries:
        if sequences_overlap(seq, targets, min_match=min_match):
            return True
    return False


def sequences_overlap(query, targets, min_match=30):
    """Check if a query overlaps with any sequence in targets."""
    for target in targets:
        cigar = align(query, target, 'HW', 'path')['cigar']
        if max([int(i) for i in re.split(r"[=IDX]+", cigar) if i]) > min_match:
            return True
        else:
            cigar = align(revcom(query), target, 'HW', 'path')['cigar']
            if cigar_to_max_operation(cigar) > min_match:
                return True
    return False


@lru_cache(maxsize=10000)
def cigar_to_max_operation(cigar):
    """Get the longest cigar operation from this cigar string."""
    return max([int(i) for i in re.split(r"[=IDX]+", cigar) if i])


@lru_cache(maxsize=10000)
def revcom(string):
    """Build reverse complement of string."""
    return "".join([COMPLEMENTARY_SEQUENCES[s] for s in string[::-1]])
