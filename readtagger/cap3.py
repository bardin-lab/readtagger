import logging
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

logger = logging.getLogger(__name__)


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
        >>> too_many_reads = {i: read1 for i in range(802)}
        >>> len(Cap3Assembly(too_many_reads).contigs)
        0
        """
        super(Cap3Assembly, self).__init__(sequences=sequences, shm_dir=shm_dir)

    def assemble(self):
        """Assemble sequences."""
        if 0 < len(self.sequences) < self.seq_limit:
            with open(os.devnull, 'w') as DEVNULL:
                args = ['cap3', self.input_path, '-p', '75', '-s', '500', '-z', '2']
                try:
                    # Use check call to ignore stdout of cap3
                    subprocess.check_call(args, stdout=DEVNULL, close_fds=True, timeout=120)
                except subprocess.SubprocessError as e:
                    logger.error("An error occured while attempting to assemble reads: "
                                 "%s\n The problematic sequences are: %s", e, self.sequences)
                    return Ace.ACEFileRecord().contigs
            return Ace.read(open(os.path.join(self.input_dir, 'multialign.fa.cap.ace'))).contigs
        else:
            # We return an empty record if there are too many sequences to assemble
            return Ace.ACEFileRecord().contigs
