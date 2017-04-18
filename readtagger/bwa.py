import os
import subprocess
import pysam
import temporary


class Bwa(object):
    """Hold Blast-related data and methods."""

    def __init__(self, input_path, bwa_index=None, reference_fasta=None, threads=1):
        """
        BWA object for `sequences`.

        Align sequences in fastq/fasta file `input_path` to bwa_index or construct a new index using reference_fasta

        >>> from test.helpers import roo_seq
        >>> from .fasta_io import write_sequences
        >>> with temporary.temp_dir() as tempdir:
        ...     reference_fasta = os.path.join(str(tempdir), 'reference.fasta')
        ...     write_sequences({'cluster_1_left_sequences_0': roo_seq}, output_path=reference_fasta)
        ...     b = Bwa(input_path=reference_fasta, reference_fasta=reference_fasta)
        >>> len(b.bwa_run) == 1
        True
        >>> result = b.reads_to_clusters()
        >>> assert len(result[1]['left_sequences']) == 1
        """
        self.input_path = input_path
        self.bwa_index = bwa_index
        self.reference_fasta = reference_fasta
        self.threads = threads
        self.bwa_run = self.run()
        self.clusters = self.reads_to_clusters()
        self.best_candidates, self.left_candidates, self.right_candidates = self.describe_clusters()
        pass

    def run(self):
        """Run bwa command."""
        with temporary.temp_dir() as temp_dir:
            temp_dir = str(temp_dir)
            if not self.bwa_index:
                self.bwa_index = make_bwa_index(self.reference_fasta, dir=temp_dir)
            proc = subprocess.Popen(['bwa', 'mem', '-t', str(self.threads), self.bwa_index, self.input_path],
                                    stdout=subprocess.PIPE, env=os.environ.copy(), close_fds=True)
            f = pysam.AlignmentFile(proc.stdout)
            self.header = f.header['SQ']
            reads = [r for r in f]
            proc.stdout.close()
            return reads

    def reads_to_clusters(self):
        """Map a readname back to a specific cluster contig or sequence."""
        clusters = {}
        for r in self.bwa_run:
            qname = r.query_name
            cluster_number = int(qname.split('cluster_')[1].split('_')[0])
            if cluster_number not in clusters:
                clusters[cluster_number] = {'left_sequences': {}, 'right_sequences': {}, 'left_contigs': {}, 'right_contigs': {}}
            for cluster_item in clusters[cluster_number].keys():
                if cluster_item in qname:
                    number = int(qname.split('_%s_' % cluster_item)[1])
                    clusters[cluster_number][cluster_item][number] = r
                    break
        return clusters

    def describe_clusters(self):
        """Return a list of possible matches, sorted by the length of the region that is being covered."""
        for cluster_number, cluster in self.clusters.items():
            left = self.split_reads_into_tid_clusters(cluster['left_contigs'] or cluster['left_sequences'])
            right = self.split_reads_into_tid_clusters(cluster['right_contigs'] or cluster['right_sequences'])
            common_tids = set(left.keys()) & set(right.keys())
            best_candidates = []
            left_candidates = []
            right_candidates = []
            if common_tids:
                for common_tid in common_tids:
                    length = self.header[common_tid]['LN']
                    min_left = min([r.pos for r in left[common_tid]])
                    min_right = min([r.pos for r in right[common_tid]])
                    max_left = max([(r.pos + r.alen) for r in left[common_tid]])
                    max_right = max([(r.pos + r.alen) for r in right[common_tid]])
                    if max_right - min_left >= max_left - min_right:
                        start = min_left
                        end = max_right
                    else:
                        start = min_right
                        end = max_left
                    full_length_fraction = (end - start) / float(length)
                    support = len(set(r.query_name for r in left[common_tid])) + len(set(r.query_name for r in right[common_tid]))
                    best_candidates.append((self.header[common_tid]['SN'], start, end, full_length_fraction, support))
            else:
                for tid, reads in left.items():
                    length = self.header[tid]['LN']
                    start = min([r.pos for r in reads])
                    end = max([(r.pos + r.alen) for r in reads])
                    full_length_fraction = (end - start) / float(length)
                    support = len(set(r.query_name for r in reads))
                    left_candidates.append((self.header[tid]['SN'], start, end, full_length_fraction, support))
                for tid, reads in right.items():
                    length = self.header[tid]['LN']
                    start = min([r.pos for r in reads])
                    end = max([(r.pos + r.alen) for r in reads])
                    full_length_fraction = (end - start) / float(length)
                    support = len(set(r.query_name for r in reads))
                    right_candidates.append((self.header[tid]['SN'], start, end, full_length_fraction, support))
            # Sort by distance between end and start. That's probably not the best idea ...
            best_candidates = sorted(best_candidates, key=lambda x: x[2] - x[1])
            left_candidates = sorted(left_candidates, key=lambda x: x[2] - x[1])
            right_candidates = sorted(right_candidates, key=lambda x: x[2] - x[1])
        return best_candidates, left_candidates, right_candidates

    def split_reads_into_tid_clusters(self, read_d):
        """Split reads in read_d into clusters based on the read tid."""
        cluster = {}
        for read in read_d.values():
            if read.tid not in cluster:
                cluster[read.tid] = [read]
            else:
                cluster[read.tid].append(read)
        return cluster


def make_bwa_index(reference_fasta, dir='.'):
    """Make a bwa index for reference_fasta and return path to the index's basename."""
    fasta_basename = os.path.basename(reference_fasta)
    target_fasta = os.path.join(dir, fasta_basename)
    os.symlink(reference_fasta, target_fasta)
    args = ['bwa', 'index', target_fasta]
    subprocess.call(args, env=os.environ.copy())
    return os.path.abspath(target_fasta)
