import os
import subprocess
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from .fasta_io import write_sequences
import temporary


class Blast(object):
    """Hold Blast-related data and methods."""

    def __init__(self, sequence, blastdb=None, reference_fasta=None):
        """
        Blast object for `sequences`.

        Writes out sequences in `sequences` dict as fasta and aligns sequences against `blastdb`.

        >>> from test.helpers import roo_seq
        >>> with temporary.temp_dir() as tempdir:
        ...     reference_fasta = os.path.join(str(tempdir), 'reference.fasta')
        ...     write_sequences({'roo_fragment': roo_seq}, output_path=reference_fasta)
        ...     s = 'ATCGAGAGCGATAAATTATATTTAGGATTTTGTTATCTAAGGCGACAGCTCAAAAACATGTAATTTAAGTGCACACTACATGAGTCAGTCACTTGAGATCGTTCCCCGCCTCCTAAAAT'
        ...     sequence = {'test': s}
        ...     b = Blast(sequence=sequence, reference_fasta=reference_fasta)
        >>> len(b.blast_result[0].alignments) == 1
        True
        """
        self.sequence = sequence
        self.blastdb = blastdb
        self.reference_fasta = reference_fasta
        self.blast_result = self.run()

    def run(self):
        """Run blast command."""
        with temporary.temp_dir() as temp_dir:
            temp_dir = str(temp_dir)
            if not self.blastdb:
                self.blastdb = make_blastdb(self.reference_fasta, dir=temp_dir)
            input_path = os.path.abspath(os.path.join(temp_dir, 'blast.fa'))
            outfile = 'blast_result.xml'
            write_sequences(sequences=self.sequence, output_path=input_path)
            cline = NcbiblastnCommandline(cmd='blastn',
                                          query=input_path,
                                          db=self.blastdb,
                                          out=outfile,
                                          evalue=0.00001,
                                          max_target_seqs=5,
                                          strand='plus',
                                          max_hsps=5,
                                          outfmt=5)
            self.stdout, self.stderr = cline()
            try:
                result = list(NCBIXML.parse(open(outfile)))
            except:
                result = []
            return result


def make_blastdb(reference_fasta, dir='.'):
    """Make a blastdb for insertion_fasta and return path to blastdb."""
    args = ['makeblastdb', '-in', reference_fasta, '-parse_seqids', '-dbtype', 'nucl', '-out', os.path.join(dir, os.path.basename(reference_fasta))]
    subprocess.call(args, env=os.environ.copy())
    return os.path.abspath(os.path.join(dir, os.path.basename(reference_fasta)))
