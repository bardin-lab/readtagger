import os
import shutil
import subprocess
import tempfile
from Bio.Sequencing import Ace


class Cap3Assembly(object):
    """A class that holds reads of a cluster and assembles them using cap3."""

    def __init__(self, sequences):
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
        self.input_dir = tempfile.mkdtemp(dir='.')
        self.input_path = os.path.join(self.input_dir, 'multialign.fa')
        self.write_sequences()
        self.assemble()
        self.assembly = Ace.read(open(os.path.join(self.input_dir, 'multialign.fa.cap.ace')))
        shutil.rmtree(self.input_dir)

    def write_sequences(self):
        """Take sequences and write them out to temporary file for cap3."""
        with open(self.input_path, 'w') as out:
            for qname, seq in self.sequences.items():
                out.write(">%s\n%s\n" % (qname, seq))

    def assemble(self):
        """Assemble sequences."""
        args = ['cap3', self.input_path, '-p', '75', '-s', '500', '-z', '2']
        subprocess.check_call(args, env=os.environ.copy())  # Use check call to ignore stdout of cap3

    @staticmethod
    def join_assemblies(assemblies):
        """Get contigs for each assembly and attempt to generate a longer contig."""
        contigs = [contig for cap3obj in assemblies for contig in cap3obj.assembly.contigs]
        sequences = {"Contig_%s" % i: contig.sequence for i, contig in enumerate(contigs)}
        return Cap3Assembly(sequences)
