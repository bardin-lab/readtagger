import pysam

IN = '/Users/mvandenb/src/tools-devteam/tools/bwa/test-data/bwa-mem-test1.bam'


class SamTagProcessor(object):
    def __init__(self, source_path, root_tag='TE'):
        self.source_alignment = pysam.AlignmentFile(source_path)
        self.result = self.process_source()

    def compute_tag(self, r):
        """
        Adds tags for downstream processing:
            - Reference
            - Pos
            - Sense
            - Aligned position
            - MD
        Returns a tag of the form:
        [('AA'), ('R:gypsy,POS:1,MD:107S35M,S:AS')]
        'AA' stands for alternative alignment
        :param r: AlignedSegement
        """
        tags = dict(ref=self.source_alignment.get_reference_name(r.tid),
                    pos=r.reference_start,
                    md=r.get_tag('MD'),
                    sense='AS' if r.is_reverse else 'S')
        return [('AA', "R:{ref},POS:{pos},MD:{md},S:{sense}".format(**tags))]

    def is_taggable(self, r):
        """
        Decide if a read should be tagged
        :param r: AlignedSegement
        :type r: pysam.Alignead
        """
        return not r.is_qcfail and not r.is_secondary and not r.is_supplementary and not r.is_unmapped

    def process_source(self):
        return {r.query_name: {'r1' if r.is_read1 else 'r2': self.compute_tag(r)} for r in self.source_alignment if self.is_taggable(r)}

    def add_mate(self):
        for qname, tag_d in self.result.items():
            if len(tag_d) == 2:
                # Both mates aligned, add mate tag
                tag_d['r1'] = tag_d['r1'].extend(('MA', tag_d['r2'][0][1]))
                tag_d['r2'] = tag_d['r2'].extend(('MA', tag_d['r1'][0][1]))
            elif len(tag_d) == 1:
                if 'r1' in tag_d:
                    tag_d['r2'] = [('MA', tag_d['r1'][0][1])]
                else:
                    tag_d['r1'] = [('MA', tag_d['r2'][0][1])]
            self.result[qname] = tag_d

class SamAnnotator(object):
    def __init__(compare_with, readtags):
        """
        Compare `readtags` with `compare_with`.
        """
        self.compare_with = pysam.AlignmentFile(source_path)
