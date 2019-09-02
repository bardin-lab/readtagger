import logging
import os

from concurrent.futures import (
    as_completed,
    wait,
    ThreadPoolExecutor,
    ProcessPoolExecutor
)

from .assemby_realignment import AssemblyRealigner
from .bam_io import (
    BamAlignmentReader as Reader,
    BamAlignmentWriter as Writer,
    merge_bam,
    sort_bam,
    split_locations_between_clusters
)
from .bwa import (
    Bwa,
    make_bwa_index
)
from .cluster import Cluster
from .cluster import collect_evidence
from .cluster_base import (
    SampleNameMixin,
    ToGffMixin,
    ToVcfMixin
)
from .cigar import aligned_segment_corresponds_to_transposable_element
from .fasta_io import merge_fasta
from .filter_insertions import sequences_match
from .find_softclip_clusters import SoftClipClusterFinder
from .gff_io import merge_gff_files
from .readtagger import get_max_proper_pair_size
from .vcf_io import merge_vcf_files
from .verify import discard_supplementary
from tempfile import (
    mkstemp,
    TemporaryDirectory,
)

logger = logging.getLogger(__name__)
DEFAULT_MAX_PROPER_PAIR_SIZE = 700


class ClusterManager(object):
    """Coordinate multiple ClusterFinder objects when running in multiprocessing mode."""

    def __init__(self, **kwds):
        """Decide if passing kwds on to ClusterFinder or if splitting input file is required."""
        if kwds['threads'] > 1:
            self.threads = kwds['threads']
            # this is ugly, but each ClusterFinder instance should be able to use an additional thread
            kwds['threads'] = 2
            self.kwds = kwds
            self.process_list = []
            self.process()
        else:
            # Delegate to clusterfinder
            ClusterFinder(**kwds)

    def process(self):
        """Process input bam in chunks."""
        with TemporaryDirectory(prefix='ClusterManager_') as tempdir:
            executor = ProcessPoolExecutor(max_workers=self.threads)
            futures = []
            chunks = split_locations_between_clusters(self.kwds['input_path'], region=self.kwds.get('region'))
            if self.kwds['transposon_reference_fasta'] and not self.kwds['transposon_bwa_index']:
                self.kwds['transposon_bwa_index'], _ = make_bwa_index(self.kwds['transposon_reference_fasta'], dir=tempdir)
            if self.kwds['genome_reference_fasta'] and not self.kwds['genome_bwa_index']:
                self.kwds['genome_bwa_index'], _ = make_bwa_index(self.kwds['genome_reference_fasta'], dir=tempdir)
            for i, region in enumerate(chunks):
                kwds = self.kwds.copy()
                kwds['region'] = region
                for key, ext in [('output_bam', '.bam'), ('output_gff', '.gff'), ('output_vcf', '.vcf'), ('output_fasta', '.fasta')]:
                    kwds[key] = os.path.join(tempdir, "%d%s" % (i, ext))
                self.process_list.append(kwds)
                futures.append(executor.submit(wrapper, kwds))
            for f in as_completed(fs=futures):
                e = f.exception()
                if e is not None:
                    if isinstance(e, RuntimeError):
                        logger.error("Runtime error occured: %s", e)
                    else:
                        logger.error("Shutting down futures, an Exception occured.")
                        wait_for_running_futures = []
                        for rf in futures[::-1]:
                            if rf.cancel():
                                wait_for_running_futures.append(rf)
                        wait(wait_for_running_futures)
                        raise f.exception()
                executor.shutdown()
            self.merge_outputs()

    def merge_outputs(self):
        """Merge outputs produced by working over smaller chunks with ClusterManager."""
        output_bam = self.kwds.get('output_bam')
        if output_bam:
            bam_files = [kwd['output_bam'] for kwd in self.process_list if kwd['output_bam']]
            merge_bam(bam_collection=bam_files, output_path=self.kwds['output_bam'], sort_order='coordinate', threads=min((8, self.threads)))
        output_fasta = self.kwds.get('output_fasta')
        if output_fasta:
            merge_fasta(fasta_files=[kwd['output_fasta'] for kwd in self.process_list], output_path=output_fasta)
        output_gff = self.kwds.get('output_gff')
        if output_gff:
            merge_gff_files([kwd['output_gff'] for kwd in self.process_list], output_gff)
        output_vcf = self.kwds.get('output_vcf')
        if output_vcf:
            merge_vcf_files([kwd['output_vcf'] for kwd in self.process_list], output_vcf)


def wrapper(kwds):
    """Launch ClusterFinder instances."""
    ClusterFinder(**kwds)


class ClusterFinder(SampleNameMixin, ToGffMixin, ToVcfMixin):
    """Find clusters of reads."""

    def __init__(self,
                 input_path,
                 output_bam=None,
                 output_gff=None,
                 output_fasta=None,
                 output_vcf=None,
                 transposon_reference_fasta=None,
                 transposon_bwa_index=None,
                 genome_reference_fasta=None,
                 genome_bwa_index=None,
                 include_duplicates=False,
                 sample_name=None,
                 threads=1,
                 min_mapq=4,
                 max_clustersupport=800,
                 max_proper_pair_size=0,
                 remove_supplementary_without_primary=False,
                 region=None,
                 shm_dir=None,
                 skip_decoy=True):
        """
        Find readclusters in input_path file.

        This class assumes that all reads in input_path that should be clustered have either an 'AD' tag or a 'BD' tag.
        The initial limits of the cluster are defined by the maximum extent of overlapping reads, and each read that is added at the 3 prime end of
        the cluster will extend the cluster.
        The join_cluster method will then join clusters that overlap through their clipped sequences and cluster that can be assembled based on their proximity
        and the fact that they support the same same insertion (and can hence contribute to the same contig if assembled).
        """
        self._sample_name = sample_name
        self.shm_dir = shm_dir
        self.region = region
        self.input_path = input_path
        self.output_bam = output_bam
        self.output_gff = output_gff
        self.output_vcf = output_vcf
        self.output_fasta = output_fasta
        self.transposon_reference_fasta = transposon_reference_fasta
        self.transposon_bwa_index = transposon_bwa_index
        self.genome_reference_fasta = genome_reference_fasta
        self.genome_bwa_index = genome_bwa_index
        self.include_duplicates = include_duplicates
        self.min_mapq = min_mapq
        self.max_clustersupport = max_clustersupport
        self.max_proper_pair_size = max_proper_pair_size or get_max_proper_pair_size(self.input_path, region=self.region) or DEFAULT_MAX_PROPER_PAIR_SIZE
        self.skip_decoy = skip_decoy
        self.softclip_finder = SoftClipClusterFinder(region=self.region,
                                                     min_mapq=self.min_mapq,
                                                     sample_name=self.sample_name)
        with TemporaryDirectory(prefix='ClusterFinder_') as self._tempdir:
            self.transposon_bwa_index, self.genome_bwa_index = self.setup_bwa_indexes()
            self.remove_supplementary_without_primary = remove_supplementary_without_primary
            self.threads = threads
            self.tp = ThreadPoolExecutor(threads)  # max threads
            if self.genome_bwa_index and self.transposon_bwa_index:
                self.assembly_realigner = AssemblyRealigner(input_alignment_file=self.input_path,
                                                            genome_bwa_index=self.genome_bwa_index,
                                                            transposon_bwa_index=self.transposon_bwa_index)
            else:
                self.assembly_realigner = None
            self.is_decoy = False
            self.clusters = self.find_cluster()
            if not self.is_decoy or not self.skip_decoy:
                self.clean_clusters()
                self.join_clusters()
                self.to_fasta()
                self.align_bwa()
                # second pass, see if having reference name assigned by align_bwa improves joining
                self.join_clusters()
                self.to_fasta()
                self.align_bwa()
                self.collect_evidence()
                self.annotate_softclip()
                self.to_bam()
                self.to_gff(output_path=self.output_gff)
                self.to_vcf(output_path=self.output_vcf)

    def setup_bwa_indexes(self):
        """Handle setting up BWA indexes."""
        transposon_bwa_index = self.transposon_bwa_index
        genome_bwa_index = self.genome_bwa_index
        if not transposon_bwa_index and self.transposon_reference_fasta:
            transposon_bwa_index, _ = make_bwa_index(self.transposon_reference_fasta, dir=self._tempdir)
        if not genome_bwa_index and self.genome_reference_fasta:
            genome_bwa_index, _ = make_bwa_index(self.genome_reference_fasta, dir=self._tempdir)
        return transposon_bwa_index, genome_bwa_index

    def _remove_supplementary_without_primary(self):
        """Remove supplementary reads without primary alignments."""
        output_path = os.path.join(self._tempdir, 'clean.bam')
        discard_supplementary(input_path=self.input_path, output_path=output_path)
        self.input_path = output_path

    def find_cluster(self):
        """Find clusters by iterating over input_path and creating clusters if reads are disjointed."""
        logger.info("Finding clusters in region '%s'", self.region or 0)
        if self.remove_supplementary_without_primary:
            self._remove_supplementary_without_primary()
        clusters = []
        skip = None
        with Reader(self.input_path, region=self.region, index=True) as reader:
            self.header = reader.header
            for r in reader.fetch(region=self.region):
                if not self.include_duplicates:
                    if r.is_duplicate:
                        continue
                if not r.mapping_quality >= self.min_mapq or r.reference_start == skip:
                    continue
                self.softclip_finder.add_read(r=r)
                if not (r.has_tag('BD') or r.has_tag('AD')):
                    continue
                if r.has_tag('AD') and not aligned_segment_corresponds_to_transposable_element(r):
                    continue
                if not clusters:
                    cluster = Cluster(shm_dir=self.shm_dir, max_proper_size=self.max_proper_pair_size)
                    cluster.append(r)
                    clusters.append(cluster)
                    continue
                if clusters[-1].read_is_compatible(r):
                    clusters[-1].append(r)
                elif len(clusters) >= 2 and r.reference_start == clusters[-1][-1].reference_start == clusters[-2][-1].reference_start:
                    skip = r.reference_start
                    clusters[-1].abnormal = True
                    clusters[-2].abnormal = True
                    # Could be X:22,432,984-22,433,240, a huge accumulation of fragments with rover homology.
                else:
                    cluster = Cluster(shm_dir=self.shm_dir, max_proper_size=self.max_proper_pair_size)
                    cluster.append(r)
                    clusters.append(cluster)
        self.softclip_finder.merge_clusters()
        logger.info('Found %d cluster on first pass (%s)', len(clusters), self.region or 0)
        if clusters:
            minimum_start = clusters[0].min
            maximum_end = clusters[-1].max
            cluster_density = len(clusters) / float(maximum_end - minimum_start)
            if cluster_density > 0.1:
                # every 10th nt a cluster, that should only happen on decoys.
                self.is_decoy = True
                logger.info('Skipping region with abnormally high cluster density (%s), probably a decoy (%s)',
                            cluster_density,
                            self.region or 0)
        return clusters

    def clean_clusters(self):
        """Remove clusters that have more reads supporting an insertion than specified in self.max_clustersupport."""
        self.clusters = [c for c in self.clusters if not len(c.read_index) > self.max_clustersupport]

    def join_clusters(self):
        """Iterate over self.cluster and attempt to join consecutive clusters."""
        new_clusterlength = 0
        if len(self.clusters) > 1:
            cluster_length = len(self.clusters)
            i = 0
            while new_clusterlength != cluster_length:
                i += 1
                logger.info("Joining clusters (currently %d), round %i (%s)", cluster_length, i, self.region or 0)
                cluster_length = new_clusterlength
                for cluster in self.clusters:
                    cluster.join_adjacent(all_clusters=self.clusters)
                new_clusterlength = len(self.clusters)
        logger.info("Found %d cluster after first pass of cluster joining (%s).", new_clusterlength, self.region or 0)
        logger.info("Splitting cluster at polarity switches")
        for index, cluster in enumerate(self.clusters):
            new_clusters = cluster.split_cluster_at_polarity_switch()
            self._add_new_clusters(new_clusters, index)
        logger.info("After splitting at polarity switches we have %d cluster (was: %d) (%s)",
                    len(self.clusters),
                    new_clusterlength,
                    self.region or 0)
        logger.info("Checking cluster consistency (%s)", self.region or 0)
        for index, cluster in enumerate(self.clusters):
            new_clusters = cluster.check_cluster_consistency()
            self._add_new_clusters(new_clusters, index)
        logger.info("After splitting inconsistent clusters we have %d cluster", len(self.clusters))
        logger.info("Last pass of joining cluster (%s)", self.region or 0)
        for cluster in self.clusters:
            cluster.refine_members(self.assembly_realigner)
            cluster.join_adjacent(all_clusters=self.clusters)
        # We are done, we can give the clusters a numeric index, so that we can distribute the processing and recover the results
        for i, cluster in enumerate(self.clusters):
            cluster.sequence = i
        logger.info("Found %d cluster overall (%s)", len(self.clusters), self.region or 0)
        self.clusters.sort(key=lambda x: x.start)

    def annotate_softclip(self):
        """Walk along all found clusters and annotate them with softclip clusters."""
        SEARCH_WINDOW = 20
        softclip_idx = 0
        for cluster in self.clusters:
            left_breakpoint_sequence = cluster.clustertag.left_breakpoint_sequence
            right_breakpoint_sequence = cluster.clustertag.right_breakpoint_sequence
            if left_breakpoint_sequence or right_breakpoint_sequence:
                start, end = cluster._start_and_end
                for i, c in enumerate(self.softclip_finder.clusters[softclip_idx:]):
                    if (start - SEARCH_WINDOW) < c.clip_position < (end + SEARCH_WINDOW):
                        if sequences_match(c.consensus, left_breakpoint_sequence, '3p_clip'):
                            cluster.softclip_clusters.append(c.id)
                        elif sequences_match(c.consensus, right_breakpoint_sequence, '5p_clip'):
                            cluster.softclip_clusters.append(c.id)
                    elif c.clip_position > end + SEARCH_WINDOW:
                        softclip_idx = i - 100
                        if softclip_idx < 0:
                            softclip_idx = 0
                        break

    def _add_new_clusters(self, new_clusters, index):
        current_index = index
        for cluster in new_clusters:
            if cluster:
                if current_index == index:
                    self.clusters[current_index] = cluster
                else:
                    self.clusters.insert(current_index, cluster)
                current_index += 1

    def collect_evidence(self):
        """Count reads that overlap cluster site but do not provide evidence for an insertion."""
        logger.info("Collecting evidence (%s)", self.region or 0)
        with Reader(self.input_path, index=True) as alignment_file:
            for cluster in self.clusters:
                if not cluster.abnormal:
                    collect_evidence(cluster, alignment_file)

    def _create_contigs(self):
        futures = []
        for cluster in self.clusters:
            futures.append(self.tp.submit(cluster._make_contigs))
        wait(futures)

    def to_fasta(self):
        """Write supporting sequences to fasta file for detailed analysis."""
        logger.info("Writing contig fasta (%s)", self.region or 0)
        if not self.output_fasta:
            fd, self.output_fasta = mkstemp(suffix='.fasta', prefix='contigs', dir=self._tempdir)
            os.close(fd)
        with open(self.output_fasta, 'w') as out:
            for cluster in self.clusters:
                for seq in cluster.to_fasta():
                    out.write(seq)

    def align_bwa(self):
        """Align cluster contigs or invidiual reads to a reference and write result into cluster."""
        logger.info("Aligning reads with BWA to describe cluster (%s)", self.region or 0)
        if self.output_fasta and (self.transposon_reference_fasta or self.transposon_bwa_index):
            bwa = Bwa(input_path=self.output_fasta,
                      bwa_index=self.transposon_bwa_index,
                      reference_fasta=self.transposon_reference_fasta,
                      threads=self.threads)
            for cluster in self.clusters:
                cluster.feature_args = []
                description = bwa.description.get(cluster.id)
                if description:
                    cluster.insert_reference_name = description.pop(-1)
                    for cluster_description in description:
                        if cluster_description:
                            cluster.feature_args.append(cluster_description[0])

    def to_bam(self):
        """Write clusters of reads and include cluster number in CD tag."""
        logger.info("Writing clusters of reads (%s)", self.region or 0)
        if self.output_bam:
            with Writer(self.output_bam, header=self.header) as writer:
                for i, cluster in enumerate(self.clusters):
                    for r in cluster:
                        r.set_tag('CD', i)
                        writer.write(r)
                    for r in cluster.evidence_for_five_p:
                        r.set_tag('CD', i)
                        r.set_tag('XD', 5)
                        writer.write(r)
                    for r in cluster.evidence_for_three_p:
                        r.set_tag('CD', i)
                        r.set_tag('XD', 3)
                        writer.write(r)
                    for r in cluster.evidence_against:
                        r.set_tag('DD', i)
                        writer.write(r)
            if self.threads < 2:
                # Because the cluster splitting doesn't necessarily conserve order we need to sort again.
                # We do this only if not using threads, since with multiple threads we split and merge at the end,
                # so there we'd need to sort again anyway.
                sort_bam(inpath=self.output_bam, output=self.output_bam, sort_order='coordinate', threads=self.threads)
