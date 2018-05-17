import re
from edlib import align

from .instance_lru import lru_cache
from .utils import revcom


def multiple_sequences_overlap(queries, targets, min_match=30, check_revcom=True):
    """Check if any query in queries overlaps with any target in targets."""
    for seq in queries:
        if sequences_overlap(seq, targets, min_match=min_match, check_revcom=check_revcom):
            return True
    return False


def sequences_overlap(query, targets, min_match=30, check_revcom=True):
    """Check if a query overlaps with any sequence in targets."""
    for target in targets:
        cigar = align(query, target, 'HW', 'path')['cigar']
        if cigar_to_max_operation(cigar) > min_match:
            return True
        elif check_revcom:
            cigar = align(revcom(query), target, 'HW', 'path')['cigar']
            if cigar_to_max_operation(cigar) > min_match:
                return True
    return False


@lru_cache(maxsize=10000)
def cigar_to_max_operation(cigar):
    """Get the longest cigar operation from this cigar string."""
    return max([int(i) for i in re.split(r"[=IDX]+", cigar) if i])
