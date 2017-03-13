import os
import shutil
import tempfile
from cached_property import cached_property
from .blast_io import (
    make_blastdb,
    Blast
)


class BlastProcessor(object):
    """Holds reference to blast database and aligns individual contigs."""

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

    def blast_cluster(self, cluster):
        """
        Blast the contigs in a cluster and produce a description of the inserted sequences.

        We want to know a few different things about the cluster that we can answer by blasting:
          - what are the best hits for the contigs ?
          - are the hits for the left and right contig compatible?
          - are there different families of TE sequences inserted?
          - how long is the subject sequences, and what fraction is covered if the insert was full length
        TODO: If there is only a contig for one side, but reads supporting the breakpoint on the other side,
        try to report also the insert at the other side based on the reads tags.
        """
        left_results = [self.blast(sequence={i: contig}) for i, contig in enumerate(cluster.left_contigs)]
        right_results = [self.blast(sequence={i: contig}) for i, contig in enumerate(cluster.right_contigs)]
        left_description = self.describe_results(results=left_results)
        right_description = self.describe_results(results=right_results)
        return Description(left_description=left_description, right_description=right_description)

    def blast(self, sequence):
        """Blast a single `sequence`."""
        return Blast(sequence=sequence, blastdb=self.blastdb)

    def describe_results(self, results):
        """Iterate over results and build a flat dictionary of possible alignments."""
        description = {}
        for blast in results:
            for blast_result in blast.blast_result:
                for alignment in blast_result.alignments:
                    description[alignment.accession] = {'length': alignment.length,
                                                        'best_expect': alignment.hsps[0].expect,
                                                        'alignments': [{'expect': hsp.expect,
                                                                        'sbjct_start': hsp.sbjct_start,
                                                                        'sbjct_end': hsp.sbjct_end} for hsp in alignment.hsps]}
        return description

    def close(self):
        """Remove tempdir."""
        self.cleanup()


class Description(object):
    """A description of a potential insert sequences."""

    def __init__(self, left_description, right_description):
        """
        Given a left_description and right_description produced by `BlastProcessor.blast_cluster` produces a description.

        Optionally a serializatiion of the blast results, trying to emphasize the best hits.
        """
        self.left_description = left_description
        self.right_description = right_description

    @cached_property
    def best_left_description(self):
        """Return best left description according to e_value."""
        return self._best_description(self.left_description)

    @cached_property
    def best_right_description(self):
        """Return best right description according to e_value."""
        return self._best_description(self.right_description)

    def _best_description(self, description):
        if not description:
            return []
        min_expect = min(v['best_expect'] for v in description.values())
        best_descriptions = []
        for reference, values in description.items():
            hsps = [alignment for alignment in values['alignments'] if alignment['expect'] == min_expect]
            hsps = [{'expect': alignment['expect'],
                     'sbjct': reference,
                     'sbjct_length': values['length'],
                     'sbjct_start': alignment['sbjct_start'],
                     'sbjct_end': alignment['sbjct_end']} for alignment in hsps]
            best_descriptions.extend(hsps)
        return best_descriptions

    @cached_property
    def best_common_references(self):
        """
        Determine most-likely inserts, based on expect value and start/end coordinate.

        To rank the inserts, we look at start and end of left alignment and start and end of right alignment among the common inserts.
        We rank the inserts by the distance between min and max of the sbjct_start and sbjct_end sites.
        We then take the sum of the e-values as the new e-values.
        """
        common_references = set(self.left_description.keys()) & set(self.right_description.keys())
        if not common_references:
            return []
        best_common_references = []
        for reference in common_references:
            left = self.left_description[reference]
            right = self.right_description[reference]
            length = left['length']
            min_left = min([alignment['sbjct_start'] for alignment in left['alignments']])
            min_left_expect = min([alignment['expect'] for alignment in left['alignments'] if alignment['sbjct_start'] == min_left])
            min_right = min([alignment['sbjct_start'] for alignment in right['alignments']])
            min_right_expect = min([alignment['expect'] for alignment in right['alignments'] if alignment['sbjct_start'] == min_right])
            max_left = max([alignment['sbjct_end'] for alignment in left['alignments']])
            max_left_expect = min([alignment['expect'] for alignment in left['alignments'] if alignment['sbjct_end'] == max_left])
            max_right = max([alignment['sbjct_end'] for alignment in right['alignments']])
            max_right_expect = min([alignment['expect'] for alignment in right['alignments'] if alignment['sbjct_end'] == max_right])
            # We've got the most extreme possible coordinates, now figure out the maximum insert that could be covered
            if max_right - min_left >= max_left - min_right:
                full_length_fraction = (max_right + min_left) / length
                combined_e = max_right_expect + min_left_expect
                start = min_left
                end = max_right
            else:
                full_length_fraction = (max_left + min_right) / length
                combined_e = max_left_expect + min_right_expect
                start = min_right
                end = max_left
            best_common_references.append({'expect': combined_e,
                                           'sbjct': reference,
                                           'sbjct_start': start,
                                           'sbjct_end': end,
                                           'sbjct_length': length,
                                           'fraction_full_length': full_length_fraction})
        min_combined_e = min([d['expect'] for d in best_common_references])
        best_common_references = [ref for ref in best_common_references if ref['expect'] == min_combined_e]
        return best_common_references

    def to_feature_args(self):
        """
        Format object for output as GFF.

        We output the reconstructed insert as well as the left and right sequences (if there are any).
        """
        feature_args = []
        for reference in self.best_common_references:
            reference['source'] = 'predicted_insertion'
            type_qual = {'type': 'predicted_insertion',
                         'qualifiers': reference}
            feature_args.append(type_qual)
        for reference in self.best_left_description:
            reference['source'] = 'left_insertion'
            type_qual = {'type': 'left_insertion',
                         'qualifiers': reference}
            feature_args.append(type_qual)
        for reference in self.best_right_description:
            reference['source'] = 'right_insertion'
            type_qual = {'type': 'right_insertion',
                         'qualifiers': reference}
            feature_args.append(type_qual)
        return feature_args
