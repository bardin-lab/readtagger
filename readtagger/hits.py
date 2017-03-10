import os
import shutil
import tempfile
from .blast_io import (
    make_blastdb,
    Blast
)


class BlastProcessor(object):
    """Holds reference to blast database and aligns individual contigs."""

    def __init__(self, blastdb=None, reference_fasta=None):
        """
        Setup a Blast class with a blast database.

        >>> from test.helpers import roo_seq
        >>> from .fasta_io import write_sequences
        >>> tempdir = tempfile.mkdtemp()
        >>> reference_fasta = os.path.join(tempdir, 'reference.fasta')
        >>> write_sequences({'roo_fragment': roo_seq}, output_path=reference_fasta)
        >>>
        >>> s = 'ATCGAGAGCGATAAATTATATTTAGGATTTTGTTATCTAAGGCGACAGCTCAAAAACATGTAATTTAAGTGCACACTACATGAGTCAGTCACTTGAGATCGTTCCCCGCCTCCTAAAAT'
        >>> sequence = {'test': s}
        >>> b = BlastProcessor(reference_fasta=reference_fasta)
        >>> result = b.blast(sequence=sequence)
        >>> len(result.blast_result[0].alignments) == 1
        True
        >>> b2 = BlastProcessor(blastdb=b.blastdb)
        >>> result = b2.blast(sequence=sequence)
        >>> len(result.blast_result[0].alignments) == 1
        True
        >>> b2.close()
        >>> b.close()
        >>> shutil.rmtree(tempdir, ignore_errors=True)
        """
        self.blastdb = blastdb
        if self.blastdb:
            self.cleanup = lambda: False
        else:
            self.cleanup = self._cleanup
            self.blastdb = make_blastdb(reference_fasta=reference_fasta, dir=tempfile.mkdtemp(dir='.'))

    def _cleanup(self):
        """Remove temporary files."""
        shutil.rmtree(os.path.split(self.blastdb)[0], ignore_errors=True)

    def blast_cluster(self, cluster):
        """
        Blast the contigs in a cluster and produce a description of the inserted sequences.

        We want to know a few different things about the cluster that we can answer by blasting:
          - what are the best hits for the contigs ?
          - are the hits for the left and right contig compatible?
          - are there different families of TE sequences inserted?
          - how long is the subject sequences, and what fraction is covered if the insert was full length
        If there is only a contig for one side, but reads supporting the breakpoint on the other side,
        try to report also the insert at the other side based on the reads tags.
        """
        left_results = [self.blast(sequence={i: contig}) for i, contig in enumerate(cluster.left_contigs)]
        right_results = [self.blast(sequence={i: contig}) for i, contig in enumerate(cluster.right_contigs)]
        left_description = self.describe_results(results=left_results)
        right_description = self.describe_results(results=right_results)
        common_inserts = self.same_insert(left_description=left_description, right_description=right_description)
        if common_inserts:
            # Figure out best supported insert with min start and max end
            for insert in common_inserts:
                pass  # TODO: implement some sort of logic to pick up the best hit given both ends

        return left_description, right_description, common_inserts

    def blast(self, sequence):
        """Blast a single `sequence`."""
        return Blast(sequence=sequence, blastdb=self.blastdb)

    def same_insert(self, left_description, right_description):
        """Check if contigs come up as aligning to the same reference."""
        left_references = set([d['reference'] for d in left_description])
        right_references = set([d['reference'] for d in right_description])
        return left_references & right_references

    def describe_results(self, results):
        description = []
        for blast in results:
            for blast_result in blast.blast_result:
                for alignment in blast_result.alignments:
                    description.append({'reference': alignment.accession,
                                        'length': alignment.length,
                                        'best_expect': alignment.hsps[0].expect,
                                        'alignments': [{'expect': hsp.expect,
                                                        'sbjct_start': hsp.sbjct_start,
                                                        'sbjct_end': hsp.sbjct_end} for hsp in alignment.hsps]})
        return description

    def close(self):
        """Remove tempdir."""
        self.cleanup()
