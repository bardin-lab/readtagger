"""Organize a cluster of interest."""
import logging
from cached_property import cached_property
from .cap3 import Cap3Assembly
from .targetsiteduplication import TargetSiteDuplication


class TagCluster(object):
    """
    Take a cluster of reads and figure out how to organize clustered reads.

    First attempts to find the Target Site Duplication (TSD). Then tries to get the breakpoint,
    and if we get the breakpoint we can extract the sequences left and right of the insertion and
    all sequences pointing into the insertion to assemble the putative inserted sequence, for the left and right side.
    """

    def __init__(self, cluster, shm_dir=None):
        """Cluster is an iterable of pysam.AlignedSegment objects."""
        self.cluster = cluster
        self.shm_dir = shm_dir
        self.tsd = TargetSiteDuplication(self.cluster)
        self.five_p_breakpoint, self.three_p_breakpoint = self.find_breakpoint()
        # TODO: Implement looking up what TE the insert(s) best match to

    @cached_property
    def left_insert(self):
        """Return insert sequence as assembled from the left side."""
        if self.left_sequences:
            if not hasattr(self, '_left_seq_cap3'):
                self._left_seq_cap3 = Cap3Assembly(self.left_sequences, shm_dir=self.shm_dir)
            return self._left_seq_cap3

    @cached_property
    def right_insert(self):
        """Return insert sequence as assembled from the right side."""
        if self.right_sequences:
            if not hasattr(self, '_right_seq_cap3'):
                self._right_seq_cap3 = Cap3Assembly(self.right_sequences, shm_dir=self.shm_dir)
            return self._right_seq_cap3

    @cached_property
    def joint_insert(self):
        """Return joint insert sequence."""
        if self.right_sequences and self.left_sequences:
            return Cap3Assembly.join_assemblies([self.left_insert, self.right_insert], shm_dir=self.shm_dir)
        return self.right_insert or self.left_insert

    def find_breakpoint(self):
        """
        Find the breakpoint of a potential insertion.

        We know that informative reads on the left of a insertion that don't overlap the insertion (i.e that are not split)
        must always be sense. Inversely, reads that are on the right of an insertion must always be antisense.
        Reads that are split should have their split region pointing into the insertion.

        Genome:                     |012345678901234567890
        Insertion:                  |----------vv----------
        Read left, Mate in TE:      |>>>>>>>>>>
        Read left, split in TE:     |    >>>>>>>>XX
        Target Site Duplication TSD:|          TT
        Read right, split in TE:    |        XX<<<<<<<<
        Read right, Mate in TE:     |           <<<<<<<<<<

        We have multiple options to find the insertion breakpoint.
         - We can identify the Target Site Duplication.
         - We can determine the closest read pairs before the sense of the aligned mates switches.
         - We can cluster the split and informative mates and then infer from the 2 (or more ...) clusters that are generated
           where the insertion should be placed.
        If we rely on TSD detection we may loose some information we can get from mate pairs.
        The last option will work for long read sequencing while potentially keeping the paired end information.
        """
        # Start by finding reads with informative splits
        if not self.tsd.is_valid:
            # Need to determine the break more closely.
            # If 3' or 5' of breakpoint has been identified -- use that.
            # TODO: Could eventually be improved by looking for softclipped positions without AD,
            # but then I would need to look at reads that are not marked to be in the cluster.
            # I could look at the BD tagged reads though and scan for softclipping ... .
            five_p, three_p = self.tsd.five_p, self.tsd.three_p
            if three_p is None:
                three_p = self.infer_three_p_from_mates()
            if five_p is None:
                five_p = self.infer_five_p_from_mates()
            if self.tsd.five_p and self.tsd.three_p and self.tsd.three_p < self.tsd.five_p:
                # An invalid TSD was found, probably because the inferred TSD is too long.
                # This could be two close-by insertions that each have support from one side,
                # and/or --less likely-- a pre-existing duplication of that sequence?
                # Is that even a single insertion? Perhaps check that evidence
                # points to the same TE? Should I split the cluster? Maybe I should just dump
                # the cluster reads as a BAM file and inspect them to see what should be done.
                warn = "Found a cluster with 5p and 3p evidence for TSD, but reads are spaced too far apart.\n"
                warn += "The cluster coordinates are tid: %s, start:%s, end: %s" % (self.cluster[0].tid, self.cluster[0].pos, self.cluster[-1].reference_end)
                logging.warn(warn)
            return five_p, three_p
        return self.tsd.five_p, self.tsd.three_p

    def infer_five_p_from_mates(self):
        """Return rightmost reference end for sense reads with BD tag."""
        five_p_reads = [r.reference_end for r in self.cluster if not r.is_reverse and r.has_tag('BD')]
        if five_p_reads:
            return max(five_p_reads)
        else:
            return None

    def infer_three_p_from_mates(self):
        """Return leftmost reference start for antisense reads with BD tag."""
        three_p_reads = [r.pos for r in self.cluster if r.is_reverse and r.has_tag('BD')]
        if three_p_reads:
            return min(three_p_reads)
        else:
            return None

    @cached_property
    def left_sequences(self):
        """
        Find reads left of a breakpoint.

        These reads need to be sense oriented if they have a BD tag, and reads with an
        AD tag should support that particular TSD end (5p for left reads).
        """
        # TODO: extend this to also return quality values
        left_sequences = {}
        for r in self.cluster:
            if r.has_tag('BD'):
                if not r.is_reverse:
                    if r.is_read1:
                        qname = "%s.1" % r.query_name
                    else:
                        qname = "%s.2" % r.query_name
                    left_sequences[qname] = r.get_tag('MS')
            if r.has_tag('AD') and r.query_name in self.tsd.five_p_support:
                if r.reference_end == self.tsd.five_p or r.pos == self.tsd.five_p:
                    left_sequences[r.query_name] = r.query_sequence
        return left_sequences

    @cached_property
    def right_sequences(self):
        """
        Find reads right of a breakpoint.

        These reads need to be antisense oriented if they have a BD tag, and reads with an
        AD tag should support that particular TSD end (3p for right reads).
        """
        # TODO: extend this to also return quality values
        right_sequences = {}
        for r in self.cluster:
            if r.has_tag('BD') and not r.has_tag('AD'):
                # Exclude reads that have both AD and BD when determining reads that support the right side,
                # because if they overlap they will already count by their AD Tag below.
                if r.is_reverse:
                    if r.is_read1:
                        qname = "%s.1" % r.query_name
                    else:
                        qname = "%s.2" % r.query_name
                    right_sequences[qname] = r.get_tag('MS')
            if r.has_tag('AD') and r.query_name in self.tsd.three_p_support:
                if r.reference_end == self.tsd.three_p or r.pos == self.tsd.three_p:
                    right_sequences[r.query_name] = r.query_sequence
        return right_sequences

    @cached_property
    def right_sequence_count(self):
        """Count all unique sequences right of a breakpoint."""
        qnames = set(qname.rsplit('.1')[0].rsplit('.2')[0] for qname in self.right_sequences)
        return len(qnames)

    @cached_property
    def left_sequence_count(self):
        """Count all unique sequences right of a breakpoint."""
        qnames = set(qname.rsplit('.1')[0].rsplit('.2')[0] for qname in self.left_sequences)
        return len(qnames)
