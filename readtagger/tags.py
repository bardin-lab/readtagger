from cached_property import cached_property
from .cigar import (
    cigartuples_to_cigarstring,
    cigar_tuple_to_cigar_length,
    cigar_to_tuple,
    cigartuples_to_named_cigartuples
)
from collections import namedtuple


class BaseTag(object):
    """Generate class template for tags."""

    def __new__(cls, header):
        """Return a Tag class that knows how to format a tag."""
        return type('NamedTagTuple', (namedtuple('tag', 'tid reference_start cigar is_reverse mapq query_alignment_start query_alignment_end'),),
                    {'__str__': lambda self: self.tag_str_template % (self.header['SQ'][self.tid]['SN'],
                                                                      self.reference_start,
                                                                      self.query_alignment_start,
                                                                      self.query_alignment_end,
                                                                      cigartuples_to_cigarstring(self.cigar),
                                                                      'AS' if self.is_reverse else 'S',
                                                                      self.mapq),
                     'reference_name': lambda self: self.header['SQ'][self.tid]['SN'],
                     'header': header,
                     'tag_str_template': "R:%s,POS:%d,QSTART:%d,QEND:%d,CIGAR:%s,S:%s,MQ:%d"})


def make_tag(template, r):
    """Return namedtuple of read attributes."""
    return template(tid=r.tid,
                    reference_start=r.reference_start,
                    cigar=r.cigar,
                    is_reverse=r.is_reverse,
                    mapq=r.mapping_quality,
                    query_alignment_start=r.query_alignment_start,
                    query_alignment_end=r.query_alignment_end)


class Tag(object):
    """Collect tag attributes and conversion."""

    def __init__(self,
                 reference_start,
                 cigar,
                 is_reverse,
                 mapq=None,
                 query_alignment_start=None,
                 query_alignment_end=None,
                 tid=None,
                 reference_name=None,
                 header=None):
        """
        Return new Tag instance from kwds.

        Note that the cigar is always wrt the reference alignment.
        When comparing Tag object by their cigar, one of the cigar needs to be inverted if the
        Tag objects are not in the same orientation.
        """
        self.reference_start = reference_start
        self._cigar = cigar
        self.is_reverse = is_reverse
        self.mapq = mapq
        self.query_alignment_start = query_alignment_start
        self.query_alignment_end = query_alignment_end
        self.tid = tid
        self._reference_name = reference_name
        self.header = header

    @cached_property
    def reference_name(self):
        """Return reference name for this instance."""
        if self._reference_name:
            return self._reference_name
        else:
            try:
                return self.header['SQ'][self.tid]['SN']
            except Exception:
                return None

    @cached_property
    def cigar_regions(self):
        """
        Return cigar regions as list of tuples in foim [(start, end), operation].

        >>> Tag(reference_start=0, cigar='20M30S', is_reverse='True', mapq=60, query_alignment_start=0, query_alignment_end=20, tid=5).cigar_regions
        [((0, 20), 0), ((20, 50), 4)]
        """
        return cigar_tuple_to_cigar_length(self.cigar)

    @cached_property
    def cigar(self):
        """
        Lazily convert cigarstring to tuple if it doesn't exist.

        >>> Tag(reference_start=0, cigar='20M30S', is_reverse='True', mapq=60, query_alignment_start=0, query_alignment_end=20, tid=5).cigar
        [CIGAR(operation=0, length=20), CIGAR(operation=4, length=30)]
        """
        if isinstance(self._cigar, str):
            self._cigar = cigar_to_tuple(self._cigar)
        self._cigar = cigartuples_to_named_cigartuples(self._cigar)
        return self._cigar

    @staticmethod
    def from_read(r, header=None):
        """
        Return Tag instance from pysam.AlignedSegment Instance.

        >>> from test.helpers import MockAlignedSegment as AlignedSegment
        >>> t = Tag.from_read(AlignedSegment(cigar='20M30S'))
        >>> isinstance(t, Tag)
        True
        >>> t.reference_name == None
        True
        """
        return Tag(tid=r.tid,
                   reference_start=r.reference_start,
                   cigar=r.cigar,
                   is_reverse=r.is_reverse,
                   mapq=r.mapping_quality,
                   query_alignment_start=r.query_alignment_start,
                   query_alignment_end=r.query_alignment_end,
                   header=header)

    @staticmethod
    def from_tag_str(tag_str):
        """
        Return Tag Instance from tag string.

        >>> t = Tag.from_tag_str('R:FBti0019061_rover_Gypsy,POS:7435,QSTART:0,QEND:34,CIGAR:34M91S,S:S,MQ:60')
        >>> isinstance(t, Tag)
        True
        >>> t.cigar == [(0, 34), (4, 91)]
        True
        >>> t = Tag.from_tag_str('R:FBti0019061_rover_Gypsy,POS:7435,QSTART:0,QEND:34,CIGAR:34M91S,S:AS,MQ:60')
        >>> t.is_reverse
        True
        """
        tag_to_attr = {'R': 'reference_name',
                       'POS': 'reference_start',
                       'QSTART': 'query_alignment_start',
                       'QEND': 'query_alignment_end',
                       'CIGAR': 'cigar',
                       'S': 'is_reverse',
                       'MQ': 'mapq'}
        integers = ['reference_start', 'query_alignment_start', 'query_alignment_end', 'mapq']
        tag_d = {tag_to_attr[k]: v for k, v in dict(item.split(':') for item in tag_str.split(',')).items()}
        if tag_d['is_reverse'] == 'S':
            tag_d['is_reverse'] = False
        else:
            tag_d['is_reverse'] = True
        for integer in integers:
            tag_d[integer] = int(tag_d.get(integer, 0))
        return Tag(**tag_d)

    def to_dict(self):
        """
        Serialize self into dictionary.

        >>> t = Tag.from_tag_str('R:FBti0019061_rover_Gypsy,POS:7435,QSTART:0,QEND:34,CIGAR:34M91S,S:S,MQ:60')
        >>> t.to_dict()['is_reverse']
        False
        """
        return {'reference_start': self.reference_start,
                'cigar': self.cigar,
                'mapq': self.mapq,
                'query_alignment_start': self.query_alignment_start,
                'query_alignment_end': self.query_alignment_end,
                'is_reverse': self.is_reverse,
                'tid': self.tid}  # Improve this by passing tid or reference name

    def to_string(self, header=None):
        """
        Serialize to tag string.

        >>> t = Tag.from_tag_str('R:FBti0019061_rover_Gypsy,POS:7435,QSTART:0,QEND:34,CIGAR:34M91S,S:S,MQ:60')
        >>> t.to_string() == 'R:FBti0019061_rover_Gypsy,POS:7435,QSTART:0,QEND:34,CIGAR:34M91S,S:S,MQ:60'
        True

        :param header:
        :return:
        """
        header = header or self.header
        return "R:%s,POS:%d,QSTART:%d,QEND:%d,CIGAR:%s,S:%s,MQ:%d" % (self.reference_name or header['SQ'][self.tid]['SN'],
                                                                      self.reference_start,
                                                                      self.query_alignment_start,
                                                                      self.query_alignment_end,
                                                                      cigartuples_to_cigarstring(self.cigar),
                                                                      'AS' if self.is_reverse else 'S',
                                                                      self.mapq)
