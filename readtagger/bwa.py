import logging
import os
import subprocess
import pysam
import temporary

from .fasta_io import write_sequences


class Bwa(object):
    """Hold Blast-related data and methods."""

    def __init__(self, input_path, bwa_index=None, reference_fasta=None, threads=1, describe_alignment=True):
        """
        BWA object for `sequences`.

        Align sequences in fastq/fasta file `input_path` to bwa_index or construct a new index using reference_fasta

        >>> from test.helpers import roo_seq
        >>> with temporary.temp_dir() as tempdir:
        ...     reference_fasta = os.path.join(str(tempdir), 'reference.fasta')
        ...     reference_fasta = write_sequences({'cluster_1_left_sequences_0': roo_seq}, output_path=reference_fasta)
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
        self.describe_alignment = describe_alignment
        self.bwa_run = self.run()
        if self.describe_alignment:
            self.clusters = self.reads_to_clusters()
            self.desciption = self.describe_clusters()

    def run(self):
        """Run bwa command."""
        with temporary.temp_dir() as temp_dir:
            temp_dir = str(temp_dir)
            if not self.bwa_index:
                self.bwa_index, _ = make_bwa_index(self.reference_fasta, dir=temp_dir)
            proc = subprocess.Popen(['bwa', 'mem', '-B9', '-O16', '-L5', '-Y', '-t', str(self.threads), self.bwa_index, self.input_path],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=os.environ.copy(), close_fds=True)
            f = pysam.AlignmentFile(proc.stdout)
            self.header = f.header
            reads = [r for r in f]
            proc.stdout.close()
            proc.wait()
            if proc.returncode != 0:
                logging.error(proc.stderr.read())
            proc.stderr.close()
            return reads

    def reads_to_clusters(self):
        """Map a readname back to a specific cluster contig or sequence."""
        clusters = {}
        for r in self.bwa_run:
            if not r.is_unmapped:
                qname = r.query_name
                cluster_number = int(qname.split('cluster_')[1].split('_')[0])
                if cluster_number not in clusters:
                    clusters[cluster_number] = {'left_sequences': {}, 'right_sequences': {}, 'left_contigs': {}, 'right_contigs': {}}
                for cluster_item in clusters[cluster_number].keys():
                    if cluster_item in qname:
                        number = qname.split('_%s_' % cluster_item)[1]
                        clusters[cluster_number][cluster_item][number] = r
                        break
        return clusters

    def describe_clusters(self):
        """Return a list of possible matches, sorted by the length of the region that is being covered."""
        cluster_description = {}
        for cluster_number, cluster in self.clusters.items():
            all_reads = {}
            all_reads['left'] = self.split_reads_into_tid_clusters(cluster['left_contigs'] or cluster['left_sequences'])
            all_reads['right'] = self.split_reads_into_tid_clusters(cluster['right_contigs'] or cluster['right_sequences'])
            common_tids = set(all_reads['left'].keys()) & set(all_reads['right'].keys())
            best_candidates = []
            left_candidates = []
            right_candidates = []
            if common_tids:
                for common_tid in common_tids:
                    length = self.header['SQ'][common_tid]['LN']
                    min_left = min([r.pos for r in all_reads['left'][common_tid]])
                    min_right = min([r.pos for r in all_reads['right'][common_tid]])
                    max_left = max([(r.pos + r.alen) for r in all_reads['left'][common_tid]])
                    max_right = max([(r.pos + r.alen) for r in all_reads['right'][common_tid]])
                    if max_right - min_left >= max_left - min_right:
                        start = min_left
                        end = max_right
                    else:
                        start = min_right
                        end = max_left
                    full_length_fraction = (end - start) / float(length)
                    support = len(set(r.query_name for r in all_reads['left'][common_tid])) + len(set(r.query_name for r in all_reads['right'][common_tid]))
                    best_candidates.append({'sbjct': self.header['SQ'][common_tid]['SN'],
                                            'sbjct_start': start,
                                            'sbjct_end': end,
                                            'fraction_full_length': full_length_fraction,
                                            'contig_support': support})
            else:
                for orientation, tid_reads in all_reads.items():
                    for tid, reads in tid_reads.items():
                        length = self.header['SQ'][tid]['LN']
                        start = min([r.pos for r in reads])
                        end = max([(r.pos + r.alen) for r in reads])
                        full_length_fraction = (end - start) / float(length)
                        support = len(set(r.query_name for r in reads))
                        if orientation == 'left':
                            left_candidates.append({'sbjct': self.header['SQ'][tid]['SN'],
                                                    'sbjct_start': start,
                                                    'sbjct_end': end,
                                                    'fraction_full_length': full_length_fraction,
                                                    'read_support': support})
                        else:
                            right_candidates.append({'sbjct': self.header['SQ'][tid]['SN'],
                                                     'sbjct_start': start,
                                                     'sbjct_end': end,
                                                     'fraction_full_length': full_length_fraction,
                                                     'read_support': support})
            # Sort by distance between end and start. That's probably not the best idea ...
            candidates = [sorted(c, key=lambda x: x['sbjct_end'] - x['sbjct_start']) for c in (best_candidates, left_candidates, right_candidates)]
            best_candidates, left_candidates, right_candidates = candidates
            reference_name = None
            if best_candidates:
                most_likely_insertion = best_candidates[0]
                reference_name = most_likely_insertion['sbjct']
                reference_name = "_".join(reference_name.split('_')[1:-1])
            else:
                if left_candidates and right_candidates:
                    left_reference_names = ["_".join(c['sbjct'].split('_')[1:-1]) for c in left_candidates]
                    right_reference_names = ["_".join(c['sbjct'].split('_')[1:-1]) for c in right_candidates]
                    overlapping_reference_names = set(left_reference_names) & set(right_reference_names)
                    if overlapping_reference_names:
                        reference_name = overlapping_reference_names.pop()  # A random overlapping name ...
                if not reference_name:
                    if left_candidates:
                        reference_name = "_".join(left_candidates[0]['sbjct'].split('_')[1:-1])
                    elif right_candidates:
                        reference_name = "_".join(right_candidates[0]['sbjct'].split('_')[1:-1])
            cluster_description[cluster_number] = self.to_feature_args(*candidates)
            cluster_description[cluster_number].append(reference_name)

        return cluster_description

    @staticmethod
    def to_feature_args(best_candidates, left_candidates, right_candidates):
        """
        Format object for output as GFF.

        We output the reconstructed insert as well as the left and right sequences (if there are any).
        """
        feature_args = []
        type_candidates = {'predicted_insertion': best_candidates, 'left_insertion': left_candidates, 'right_insertion': right_candidates}
        for type, candidates in type_candidates.items():
            for reference in candidates:
                reference['source'] = type
                type_qual = {'type': type,
                             'qualifiers': reference}
                feature_args.append(type_qual)
        return feature_args

    def split_reads_into_tid_clusters(self, read_d):
        """Split reads in read_d into clusters based on the read tid."""
        cluster = {}
        for read in read_d.values():
            if read.tid not in cluster:
                cluster[read.tid] = [read]
            else:
                cluster[read.tid].append(read)
        return cluster


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


def make_bwa_index(reference_fasta, dir='.'):
    """Make a bwa index for reference_fasta and return path to the index's basename."""
    fasta_basename = os.path.basename(reference_fasta)
    target_fasta = os.path.abspath(os.path.join(dir, fasta_basename))
    if not os.path.exists(target_fasta):
        os.symlink(os.path.abspath(reference_fasta), target_fasta)
    args = ['bwa', 'index', target_fasta]
    p = subprocess.Popen(args, env=os.environ.copy(), close_fds=True, stderr=subprocess.PIPE)
    p.wait()
    if p.returncode != 0:
        logging.warning(p.stderr.read())
    return target_fasta, p.returncode
