import copy
import os
import shutil
import subprocess
import tempfile
from Bio.Sequencing import Ace
from .fasta_io import write_sequences


class Cap3Assembly(object):
    """A class that holds reads of a cluster and assembles them using cap3."""

    def __init__(self, sequences, shm_dir=None):
        """Asssemble sequences into contigs.

        :param sequences: dictionary with query_name as key and read sequence as value
        :type sequences: dictionary

        >>> read1 = 'TAGTTGTAAGCGATTCTTAACTTACCTACCTACATATATATACTTACGTATTTTACTATT'
        >>> read2 = 'CGAGTCGAACAAATGATCCGTCGTTTGACTAAGATCAACGCCTTTAAAGAAGTTTCAGAA'
        >>> read3 = 'TACCTACCTACATATATATACTTACGTATTTTACTATTCGAGTCGAACAAATGATCCGTC'
        >>> read4 = 'CGATTCTTAACTTACCTACCTACATATATATACTTACGTATTTTACTATTCGAGTCGAACA'
        >>> sequences = {'read1': read1, 'read2': read2, 'read3': read3, 'read4': read4}
        >>> Cap3Assembly(sequences).assembly.ncontigs
        1
        """
        self.sequences = sequences
        self.input_dir = tempfile.mkdtemp(dir=shm_dir)
        self.input_path = os.path.join(self.input_dir, 'multialign.fa')
        self.write_sequences()
        self.assemble()
        self.assembly = Ace.read(open(os.path.join(self.input_dir, 'multialign.fa.cap.ace')))
        shutil.rmtree(self.input_dir, ignore_errors=True)

    def write_sequences(self):
        """Take sequences and write them out to a temporary file for cap3."""
        write_sequences(sequences=self.sequences, output_path=self.input_path)

    def assemble(self):
        """Assemble sequences."""
        with open(os.devnull, 'w') as DEVNULL:
            args = ['cap3', self.input_path, '-p', '75', '-s', '500', '-z', '2']
            subprocess.check_call(args, stdout=DEVNULL, env=os.environ.copy(), close_fds=True)  # Use check call to ignore stdout of cap3

    @staticmethod
    def join_assemblies(assemblies, shm_dir=None):
        """Get contigs for each assembly and attempt to generate a longer contig."""
        contigs = [contig for cap3obj in assemblies for contig in cap3obj.assembly.contigs]
        sequences = {"Contig_%s" % i: contig.sequence for i, contig in enumerate(contigs)}
        return Cap3Assembly(sequences, shm_dir=shm_dir)

    @staticmethod
    def sequences_contribute_to_same_contig(seq1, seq2):
        """
        Check if a contig can be generated from two sets of sequences where both sets of sequences contributed to the contig.

        >>> read1 = 'TAGTTGTAAGCGATTCTTAACTTACCTACCTACATATATATACTTACGTATTTTACTATT'
        >>> read2 = 'CGAGTCGAACAAATGATCCGTCGTTTGACTAAGATCAACGCCTTTAAAGAAGTTTCAGAA'
        >>> read3 = 'TACCTACCTACATATATATACTTACGTATTTTACTATTCGAGTCGAACAAATGATCCGTC'
        >>> read4 = 'CGATTCTTAACTTACCTACCTACATATATATACTTACGTATTTTACTATTCGAGTCGAACA'
        >>> read5 = 'ATGCATGCATGCATGCATGCATGCATGCATGCAAAAAA'
        >>> read6 = 'ATGCATGCATGCATGCATGCATGCATGCATGCTTTTTT'
        >>> Cap3Assembly.sequences_contribute_to_same_contig(seq1={'seq1': read1, 'seq2': read2}, seq2={'seq3': read3, 'seq4': read4})
        True
        >>> Cap3Assembly.sequences_contribute_to_same_contig(seq1={'seq1': read1, 'seq2': read2}, seq2={'seq5': read5, 'seq6': read6})
        False
        """
        sequences = copy.deepcopy(seq1)
        sequences.update(seq2)
        seq1_keys = set(seq1.keys())
        seq2_keys = set(seq2.keys())
        assembly = Cap3Assembly(sequences=sequences)
        for contig in assembly.assembly.contigs:
            contig_reads = set([read.rd.name for read in contig.reads])
            if seq1_keys & contig_reads and seq2_keys & contig_reads:
                # both right_reads and left_reads contribute to the same contig,
                # this cluster should clearly be merged.
                return True
        return False
