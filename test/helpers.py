from collections import namedtuple
from itertools import chain

MockAlignedSegmentTemplate = namedtuple('AlignedSegment',
                                        ['cigar',
                                         'query_name',
                                         'reference_start',
                                         'reference_end',
                                         'pos',
                                         'has_tag',
                                         'query_alignment_start',
                                         'query_alignment_end',
                                         'query_length',
                                         'tid',
                                         'is_reverse',
                                         'mapping_quality'])


def namedtuple_to_argv(nt, prog='prog.py'):
    """Convert a namedtuple into an argv list for testing argparser arguments."""
    d = nt._asdict()
    keys = ["--%s" % k for k in d]
    argv = list(chain.from_iterable(zip(keys, d.values())))
    argv = [arg for arg in argv if not isinstance(arg, bool)]
    argv = [arg if not isinstance(arg, list) else " ".join(arg) for arg in argv]
    argv.insert(0, prog)
    return argv


class MockAlignedSegment(MockAlignedSegmentTemplate):
    """A class that simulated pysam's AlignedSegment class."""

    def __new__(cls, cigar='125M',
                query_name='read1',
                reference_start=0,
                reference_end=125,
                pos=0,
                has_tag=lambda x: True,
                query_alignment_start=0,
                query_alignment_end=125,
                query_length=125,
                tid=0,
                is_reverse=False,
                mapping_quality=60):
        """Initiate a new AlignedSegment."""
        # add default values
        return super(MockAlignedSegment, cls).__new__(cls,
                                                      cigar,
                                                      query_name,
                                                      reference_start,
                                                      reference_end,
                                                      pos,
                                                      has_tag,
                                                      query_alignment_start,
                                                      query_alignment_end,
                                                      query_length,
                                                      tid,
                                                      is_reverse,
                                                      mapping_quality)
