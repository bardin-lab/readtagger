import pysam

from .bam_io import index_bam
from .bwa import SimpleAligner
from .cap3 import Cap3Assembly
from .cigar import alternative_alignment_cigar_is_better
from .readtagger import SamTagProcessor


class AssemblyRealigner(object):
    """Assemble and realign reads in a cluster."""

    def __init__(self,
                 input_alignment_file,
                 genome_bwa_index,
                 transposon_bwa_index,
                 tmp_dir=None):  # Make that configurable ...
        """Assemble reads and align contigs in a cluster to improve breakpoints."""
        self.input_alignment_file = input_alignment_file
        self.reference_genome_index = genome_bwa_index
        self.transposon_index = transposon_bwa_index
        self.genome_aligner = SimpleAligner(bwa_index=genome_bwa_index, tmp_dir=tmp_dir)
        self.transposon_aligner = SimpleAligner(bwa_index=transposon_bwa_index, tmp_dir=tmp_dir)

    def collect_reads(self, cluster):
        """Collect reads that could be useful for determining a clusters breakpoint via assembly."""
        if len(cluster.orientation_switches) > 1:
            # For now only use this strategy to refine potential TSDs,
            # may be interesting for improving insertions with single-sided evidence as well.
            start = cluster.min
            end = cluster.max
        else:
            return []
        index_bam(self.input_alignment_file)  # Should not be necessary, but tests fail without this :(
        additional_read_gen = pysam.AlignmentFile(self.input_alignment_file).fetch(tid=cluster.tid, start=start, end=end)
        additional_reads = (r for r in additional_read_gen if r.query_name not in cluster.read_index)
        reads = {}
        for r in additional_reads:
            if r.has_tag('MS') and not r.is_duplicate:
                qname = r.query_name
                if r.is_read1:
                    current_qname = "%s.1" % qname
                    other_qname = "%s.2" % qname
                else:
                    current_qname = "%s.2" % qname
                    other_qname = "%s.1" % qname
                reads[current_qname] = r.query_sequence
                reads[other_qname] = r.get_tag('MS')
        return self.assemble_reads(reads)

    def assemble_reads(self, reads):
        """Assemble potentially informative reads, align and set tags for contigs."""
        if len(reads) > 500:
            return []
        assembly = Cap3Assembly(reads)
        contig_sequences = {i: contig.sequence for i, contig in enumerate(assembly.assembly.contigs)}
        genome_aligned_contigs, genome_header = self.genome_aligner.align_contigs(contig_sequences)
        transposon_aligned_contigs, transposon_header = self.transposon_aligner.align_contigs(contig_sequences)
        informative_reads = []
        if transposon_aligned_contigs and genome_aligned_contigs:
            transposon_tags = SamTagProcessor(source_bam=transposon_aligned_contigs, header=transposon_header, tag_mate=False)
            for gc in genome_aligned_contigs:
                if gc.query_name in transposon_tags.result:
                    tag = transposon_tags.result[gc.query_name][False]['s']
                    if alternative_alignment_cigar_is_better(current_cigar=gc.cigar,
                                                             alternative_cigar=tag.cigar,
                                                             same_orientation=gc.is_reverse == tag.is_reverse):
                        gc.set_tag('AD', str(tag))
                        gc.set_tag('AR', str(tag.reference_name()))
                        gc.set_tag('AC', gc.query_name)
                        informative_reads.append(gc)
        return informative_reads
