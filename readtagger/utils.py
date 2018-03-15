from .instance_lru import lru_cache

COMPLEMENTARY_SEQUENCES = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}


@lru_cache(maxsize=10000)
def revcom(string):
    """Build reverse complement of string."""
    return "".join([COMPLEMENTARY_SEQUENCES[s] for s in string[::-1]])
