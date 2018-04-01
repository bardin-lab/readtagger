import logging
import os
import subprocess
import pysam
from collections import defaultdict
try:
    from tempfile import TemporaryDirectory
except ImportError:
    from backports.tempfile import TemporaryDirectory

from .fasta_io import write_sequences

logger = logging.getLogger(__name__)


class Description(object):
    """Hold details for alignment."""

    def __init__(self,
                 sbjct,
                 sbjct_start,
                 sbjct_end,
                 type=None,
                 fraction_full_length=None,
                 contig_support=None,
                 read_support=None):
        """Hold and serialize description of an insertion."""
        self.sbjct = sbjct
        self.sbjct_start = sbjct_start
        self.sbjct_end = sbjct_end
        self.type = type
        self.fraction_full_length = fraction_full_length
        self.contig_support = contig_support
        self.read_support = read_support

    def to_feature_args(self):
        """Return this descriptions attributes as a GFF subfeature."""
        return {'qualifiers': {attr: value for attr, value in self.__dict__.items() if value},
                'type': self.type}


class Bwa(object):
    """Hold alignment data and methods."""

    def __init__(self, input_path, bwa_index=None, reference_fasta=None, threads=1, describe_alignment=True):
        """
        BWA object for `sequences`.

        Align sequences in fastq/fasta file `input_path` to bwa_index or construct a new index using reference_fasta

        >>> from tests.helpers import roo_seq
        >>> with TemporaryDirectory(prefix='bwa_doctest') as tempdir:
        ...     reference_fasta = os.path.join(str(tempdir), 'reference.fasta')
        ...     reference_fasta = write_sequences({'1|lsequences|0': roo_seq}, output_path=reference_fasta)
        ...     b = Bwa(input_path=reference_fasta, reference_fasta=reference_fasta)
        >>> len(b.bwa_run) == 1
        True
        >>> result = b.reads_to_clusters()
        >>> assert len(result['1']['lsequences']) == 1
        >>> b.cleanup_index()
        """
        self.input_path = input_path
        self.bwa_index = bwa_index
        self.reference_fasta = reference_fasta
        self.threads = threads
        self.describe_alignment = describe_alignment
        self.bwa_run = self.run()
        if self.describe_alignment:
            self.clusters = self.reads_to_clusters()
            self.description = self.describe_clusters()

    def run(self):
        """Run bwa command."""
        with TemporaryDirectory(prefix='BWA') as temp_dir:
            temp_dir = str(temp_dir)
            if not self.bwa_index:
                self.bwa_index, _ = make_bwa_index(self.reference_fasta, dir=temp_dir)
            proc = subprocess.Popen(['bwa', 'mem', '-B9', '-O16', '-L5', '-Y', '-t', str(self.threads), self.bwa_index, self.input_path],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
            f = pysam.AlignmentFile(proc.stdout)
            self.header = f.header
            reads = [r for r in f]
            proc.stdout.close()
            wait_and_get_return_code(proc)
            proc.stderr.close()
            return reads

    def reads_to_clusters(self):
        """Map a readname back to a specific cluster contig or sequence."""
        clusters = {}
        for r in self.bwa_run:
            if not r.is_unmapped:
                qname = r.query_name
                cluster_number = qname.split('|')[0]
                if cluster_number not in clusters:
                    clusters[cluster_number] = {'lsequences': {}, 'rsequences': {}, 'lcontigs': {}, 'rcontigs': {}}
                for cluster_item in clusters[cluster_number].keys():
                    if cluster_item in qname:
                        number = qname.split('|')[0]
                        clusters[cluster_number][cluster_item][number] = r
                        break
        return clusters

    def describe_clusters(self):
        """Return a list of possible matches, sorted by the length of the region that is being covered."""
        cluster_description = {}
        for cluster_number, cluster in self.clusters.items():
            all_reads = {}
            all_reads['left'] = self.split_reads_into_tid_clusters(cluster['lcontigs'] or cluster['lsequences'])
            all_reads['right'] = self.split_reads_into_tid_clusters(cluster['rcontigs'] or cluster['rsequences'])
            common_tids = set(all_reads['left'].keys()) & set(all_reads['right'].keys())
            best_candidates = []
            left_candidates = []
            right_candidates = []
            if common_tids:
                for common_tid in common_tids:
                    length = self.header['SQ'][common_tid]['LN']
                    min_left = min([r.reference_start for r in all_reads['left'][common_tid]])
                    min_right = min([r.reference_start for r in all_reads['right'][common_tid]])
                    max_left = max([(r.reference_start + r.reference_length) for r in all_reads['left'][common_tid]])
                    max_right = max([(r.reference_start + r.reference_length) for r in all_reads['right'][common_tid]])
                    if max_right - min_left >= max_left - min_right:
                        start = min_left
                        end = max_right
                    else:
                        start = min_right
                        end = max_left
                    full_length_fraction = (end - start) / float(length)
                    support = len(set(r.query_name for r in all_reads['left'][common_tid])) + len(set(r.query_name for r in all_reads['right'][common_tid]))
                    best_candidates.append(Description(sbjct=self.header.references[common_tid],
                                                       sbjct_start=start,
                                                       sbjct_end=end,
                                                       type='predicted_insertion',
                                                       fraction_full_length=full_length_fraction,
                                                       contig_support=support))
            else:
                for orientation, tid_reads in all_reads.items():
                    for tid, reads in tid_reads.items():
                        length = self.header.lengths[tid]
                        start = min([r.reference_start for r in reads])
                        end = max([(r.reference_start + r.reference_length) for r in reads])
                        full_length_fraction = (end - start) / float(length)
                        support = len(set(r.query_name for r in reads))
                        candidates = left_candidates if orientation == 'left' else right_candidates
                        candidates.append(Description(sbjct=self.header.references[tid],
                                                      sbjct_start=start,
                                                      sbjct_end=end,
                                                      type="%s_insert" % orientation,
                                                      fraction_full_length=full_length_fraction,
                                                      read_support=support))
            # Sort by distance between end and start. That's probably not the best idea ...
            candidates = [sorted(c, key=lambda x: -(x.sbjct_end - x.sbjct_start)) for c in (best_candidates, left_candidates, right_candidates)]
            best_candidates, left_candidates, right_candidates = candidates
            reference_name = None
            if best_candidates:
                most_likely_insertion = best_candidates[0]
                reference_name = split_reference_name(most_likely_insertion.sbjct)
            else:
                if left_candidates and right_candidates:
                    left_reference_names = [split_reference_name(c.sbjct) for c in left_candidates]
                    right_reference_names = [split_reference_name(c.sbjct) for c in right_candidates]
                    overlapping_reference_names = set(left_reference_names) & set(right_reference_names)
                    if overlapping_reference_names:
                        reference_name = overlapping_reference_names.pop()  # A random overlapping name ...
                if not reference_name:
                    candidates = left_candidates or right_candidates
                    if candidates:
                        reference_name = split_reference_name(candidates[0].sbjct)
            cluster_description[cluster_number] = [best_candidates, left_candidates, right_candidates, reference_name]
        return cluster_description

    def split_reads_into_tid_clusters(self, read_d):
        """Split reads in read_d into clusters based on the read tid."""
        cluster = defaultdict(list)
        for read in read_d.values():
            if read.tid not in cluster:
                cluster[read.tid].append(read)
        return cluster

    def cleanup_index(self):
        """Remove this instances index files."""
        cleanup_index(self.bwa_index)


class SimpleAligner(object):
    """Perform simple alignments, e.g to see if a read is contained in a contig."""

    def __init__(self, reference_sequences=None, bwa_index=None, tmp_dir=None):
        """Perform simple alignments, e.g to see if a read is contained in a contig."""
        self.tmp_dir = tmp_dir
        if not bwa_index:
            self.reference_fasta = write_sequences(reference_sequences)
            self.index, self.return_code = make_bwa_index(self.reference_fasta)
        else:
            self.index = bwa_index

    def __enter__(self):  # noqa: D105
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Clean up index when class is used as context manager."""
        self.cleanup_index()

    def align(self, sequence):
        """Return contig numbers with valid alignments for sequence."""
        sequences = write_sequences(sequence, tmp_dir=self.tmp_dir)
        aligned_reads = Bwa(input_path=sequences, bwa_index=self.index, describe_alignment=False)
        return set(r.tid for r in aligned_reads.bwa_run if not r.is_unmapped)

    def align_contigs(self, contigs):
        """Align contigs and return aligned results and header."""
        sequences = write_sequences(contigs, tmp_dir=self.tmp_dir)
        bwa_result = Bwa(input_path=sequences, bwa_index=self.index, describe_alignment=False)
        return bwa_result.bwa_run, bwa_result.header

    def cleanup_index(self):
        """Remove this instances index files."""
        cleanup_index(self.index)


def cleanup_index(index_path):
    """Remove indexes."""
    index_suffixes = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    index_files = ["%s%s" % (index_path, sfx) for sfx in index_suffixes]
    for idx_file in index_files + [index_path]:
        try:
            os.remove(idx_file)
        except Exception:
            pass


def make_bwa_index(reference_fasta, dir='.'):
    """Make a bwa index for reference_fasta and return path to the index's basename."""
    fasta_basename = os.path.basename(reference_fasta)
    target_fasta = os.path.abspath(os.path.join(dir, fasta_basename))
    if not os.path.exists(target_fasta):
        os.symlink(os.path.abspath(reference_fasta), target_fasta)
    args = ['bwa', 'index', target_fasta]
    return_code = wait_and_get_return_code(subprocess.Popen(args,
                                                            env=os.environ.copy(),
                                                            close_fds=True,
                                                            stderr=subprocess.PIPE))
    return target_fasta, return_code


def wait_and_get_return_code(p):
    """
    Wait for subprocess p, log warning on non-zero rexit code and return exit code.

    >>> p = subprocess.Popen('echo bla 1>&2 && exit 1', stderr=subprocess.PIPE, shell=True)
    >>> assert wait_and_get_return_code(p) == 1
    """
    rc = p.wait()
    if rc:
        logger.warning(p.stderr.read())
    return rc


def split_reference_name(reference_name, at="_"):
    """
    Split reference_name at the first `at`.

    This returns the TE family name if the TEs are named "$ID_$NAME_$SUPERFAMILY".

    >>> split_reference_name('FBti0019298_rover_Gypsy')
    'rover'
    """
    if at in reference_name:
        return "_".join(reference_name.split(at)[1:-1])
    else:
        return reference_name
