import os
import shutil
import subprocess
import tempfile
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from .fasta_io import write_sequences


class Blast(object):
    """Hold Blast-related data and methods."""

    def __init__(self, sequences, blastdb):
        """
        Blast object for `sequences`.

        Writes out sequences in `sequences` dict as fasta and aligns sequences against `blastdb`.

        >>> from test.helpers import roo_seq
        >>> tempdir = tempfile.mkdtemp()
        >>> reference_fasta = os.path.join(tempdir, 'reference.fasta')
        >>> write_sequences({'roo_fragment': roo_seq}, output_path=reference_fasta)
        >>> db = make_blastdb(reference_fasta, dir=tempdir)
        >>> s = 'ATCGAGAGCGATAAATTATATTTAGGATTTTGTTATCTAAGGCGACAGCTCAAAAACATGTAATTTAAGTGCACACTACATGAGTCAGTCACTTGAGATCGTTCCCCGCCTCCTAAAAT'
        >>> sequences = {'test': s}
        >>> b = Blast(sequences=sequences, blastdb=db)
        >>> len(b.blast_result[0].alignments) == 1
        True
        >>> shutil.rmtree(tempdir, ignore_errors=True)
        """
        self.sequences = sequences
        self.blastdb = blastdb
        self.input_dir = tempfile.mkdtemp(dir='.')
        self.input_path = os.path.join(self.input_dir, 'blast.fa')
        self.blast_result = self.run()
        shutil.rmtree(self.input_dir, ignore_errors=True)

    def run(self):
        """Run blast command."""
        outfile = os.path.join(self.input_dir, 'blast_result.xml')
        write_sequences(sequences=self.sequences, output_path=self.input_path)
        cline = NcbiblastnCommandline(cmd='blastn',
                                      query=self.input_path,
                                      db=self.blastdb,
                                      out=outfile,
                                      evalue=0.00001,
                                      max_target_seqs=2,
                                      max_hsps=2,
                                      outfmt=5)
        self.stdout, self.stderr = cline()
        return list(NCBIXML.parse(open(outfile)))


def make_blastdb(reference_fasta, dir='.'):
    """Make a blastdb for insertion_fasta and return path to blastdb."""
    args = ['makeblastdb', '-in', reference_fasta, '-parse_seqids', '-dbtype', 'nucl', '-out', os.path.join(dir, reference_fasta)]
    subprocess.call(args, env=os.environ.copy())
    return reference_fasta
