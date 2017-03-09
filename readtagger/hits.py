import os
import shutil
import tempfile
from .blast_io import (
    make_blastdb,
    Blast
)


class BlastProcessor(object):
    """Holds reference to blast database and aligns individual sequences."""

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

    def blast(self, sequence):
        """Blast a single `sequence`."""
        return Blast(sequence=sequence, blastdb=self.blastdb)

    def close(self):
        """Remove tempdir."""
        self.cleanup()
