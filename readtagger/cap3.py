import os
import subprocess
import abc

from Bio.Sequencing import Ace
from .fasta_io import write_sequences

# compatible with Python 2 *and* 3:
ABC = abc.ABCMeta('ABC', (object,), {'__slots__': ()})
try:
    from tempfile import TemporaryDirectory
except ImportError:
    from backports.tempfile import TemporaryDirectory


class BaseAssembly(ABC):
    """Provide Base Class for Assembly modules."""

    def __init__(self, sequences, shm_dir):
        """Run assembly."""
        self.sequences = sequences
        with TemporaryDirectory(prefix="%s" % type(self).__name__, dir=shm_dir) as self.input_dir:
            self.input_path = os.path.join(self.input_dir, 'multialign.fa')
            self.write_sequences()
            self.contigs = self.assemble()

    @abc.abstractmethod
    def assemble(self):
        """Must return contigs."""
        pass

    def write_sequences(self):
        """Take sequences and write them out to a temporary file for cap3."""
        write_sequences(sequences=self.sequences, output_path=self.input_path)


class Cap3Assembly(BaseAssembly):
    """A class that holds reads of a cluster and assembles them using cap3."""

    seq_limit = 800

    def __init__(self, sequences, shm_dir=None):
        """Asssemble sequences into contigs.

        :param sequences: dictionary with query_name as key and read sequence as value
        :type sequences: dictionary

        >>> read1 = 'TAGTTGTAAGCGATTCTTAACTTACCTACCTACATATATATACTTACGTATTTTACTATT'
        >>> read2 = 'CGAGTCGAACAAATGATCCGTCGTTTGACTAAGATCAACGCCTTTAAAGAAGTTTCAGAA'
        >>> read3 = 'TACCTACCTACATATATATACTTACGTATTTTACTATTCGAGTCGAACAAATGATCCGTC'
        >>> read4 = 'CGATTCTTAACTTACCTACCTACATATATATACTTACGTATTTTACTATTCGAGTCGAACA'
        >>> sequences = {'read1': read1, 'read2': read2, 'read3': read3, 'read4': read4}
        >>> len(Cap3Assembly(sequences).contigs)
        1
        """
        super(Cap3Assembly, self).__init__(sequences=sequences, shm_dir=shm_dir)

    def assemble(self):
        """Assemble sequences."""
        if len(self.sequences) < self.seq_limit:
            with open(os.devnull, 'w') as DEVNULL:
                args = ['cap3', self.input_path, '-p', '75', '-s', '500', '-z', '2']
                subprocess.check_call(args, stdout=DEVNULL, env=os.environ.copy(), close_fds=True)  # Use check call to ignore stdout of cap3
            return Ace.read(open(os.path.join(self.input_dir, 'multialign.fa.cap.ace'))).contigs
        else:
            # We return an empty record if there are too many sequences to assemble
            return Ace.ACEFileRecord()

    @staticmethod
    def join_assemblies(assemblies, shm_dir=None):
        """Get contigs for each assembly and attempt to generate a longer contig."""
        contigs = [contig for cap3obj in assemblies for contig in cap3obj.contigs]
        sequences = {"Contig_%s" % i: contig.sequence for i, contig in enumerate(contigs)}
        return Cap3Assembly(sequences, shm_dir=shm_dir)
