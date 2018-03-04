import logging
import os
import shutil

from cached_property import cached_property
from concurrent.futures import (
    as_completed,
    wait,
    ThreadPoolExecutor,
    ProcessPoolExecutor
)
import multiprocessing_logging
import six

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
from .cluster import non_evidence
from .gff_io import (
    sort_gff,
    write_cluster
)
from .readtagger import get_max_proper_pair_size
from .verify import discard_supplementary
try:
    from tempfile import TemporaryDirectory
except ImportError:
    from backports.tempfile import TemporaryDirectory

logger = logging.getLogger(__name__)
logging.basicConfig(format='%(asctime)s %(name)s %(levelname)s - %(message)s', level=logging.DEBUG)
multiprocessing_logging.install_mp_handler()


class ClusterManager(object):
    """Coordinate multiple ClusterFinder objects when running in multiprocessing mode."""

    def __init__(self, **kwds):
        """Decide if passing kwds on to ClusterFinder or if splitting input file is required."""
        if kwds.get('max_proper_pair_size', 0) == 0:
            kwds['max_proper_pair_size'] = get_max_proper_pair_size(kwds['input_path'])
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
                kwds['output_bam'] = os.path.join(tempdir, "%d.bam" % i)
                kwds['output_gff'] = os.path.join(tempdir, "%d.gff" % i)
                kwds['output_fasta'] = os.path.join(tempdir, "%d.fasta" % i)
                self.process_list.append(kwds)
                futures.append(executor.submit(wrapper, kwds))
            for f in as_completed(fs=futures):
                e = f.exception()
                if e is not None:
                    if isinstance(e, RuntimeError):
                        logging.error("Runtime error occured: %s", e)
                    else:
                        logging.error("Shutting down futures, an Exception occured.")
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
            merge_bam(bam_collection=bam_files, output_path=self.kwds['output_bam'])
            sort_bam(inpath=self.kwds['output_bam'], output=self.kwds['output_bam'], sort_order='coordinate', threads=min((8, self.threads)))
        output_fasta = self.kwds.get('output_fasta')
        if output_fasta:
            fasta_files = [kwd['output_fasta'] for kwd in self.process_list if kwd['output_fasta']]
            with open(output_fasta, 'wb') as fasta_writer:
                for fasta in fasta_files:
                    if os.path.exists(fasta):
                        shutil.copyfileobj(open(fasta, 'rb'), fasta_writer)
        output_gff = self.kwds.get('output_gff')
        if output_gff:
            gff_files = [kwd['output_gff'] for kwd in self.process_list if kwd['output_gff']]
            wrote_header = False
            with open(output_gff, 'w') as gff_writer:
                for gff in gff_files:
                    if os.path.exists(gff):
                        with open(gff) as gff_file:
                            for line in gff_file:
                                if line.startswith('#') and wrote_header:
                                    continue
                                gff_writer.write(line)
                        wrote_header = True
            sort_gff(input_path=output_gff, output_path=output_gff)


def wrapper(kwds):
    """Launch ClusterFinder instances."""
    ClusterFinder(**kwds)


class SampleNameMixin(object):
    """Provide a sample name property."""

    @cached_property
    def sample_name(self):
        """Return sample name if passed in manually, else guess sample name from input file."""
        if not self._sample_name:
            basename = os.path.basename(self.input_path)
            if '.' in basename:
                basename = basename.rsplit('.', 1)[0]
            return basename
        else:
            return self._sample_name


class ToGffMixin(object):
    """Provide a `to_gff` function."""

    def to_gff(self):
        """Write clusters as GFF file."""
        logging.info("Writing clusters of GFF (%s)", self.region or 0)
        if self.output_gff:
            write_cluster(clusters=self.clusters,
                          header=self.header,
                          output_path=self.output_gff,
                          sample=self.sample_name,
                          threads=self.threads)
            if self.threads < 2:
                sort_gff(self.output_gff, output_path=self.output_gff)


class ClusterFinder(SampleNameMixin, ToGffMixin):
    """Find clusters of reads."""

    def __init__(self,
                 input_path,
                 output_bam=None,
                 output_gff=None,
                 output_fasta=None,
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
        self.output_fasta = output_fasta
        self.transposon_reference_fasta = transposon_reference_fasta
        self.transposon_bwa_index = transposon_bwa_index
        self.genome_reference_fasta = genome_reference_fasta
        self.genome_bwa_index = genome_bwa_index
        self.include_duplicates = include_duplicates
        self.min_mapq = min_mapq
        self.max_clustersupport = max_clustersupport
        self.max_proper_pair_size = max_proper_pair_size
        self.skip_decoy = skip_decoy
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
                self.collect_non_evidence()
                self.to_bam()
                self.to_gff()

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
        logging.info("Finding clusters in region '%s'", self.region or 0)
        if self.remove_supplementary_without_primary:
            self._remove_supplementary_without_primary()
        clusters = []
        with Reader(self.input_path, external_bin=False, region=self.region, index=True) as reader:
            self.header = reader.header
            for r in reader.fetch(region=self.region):
                if not self.include_duplicates:
                    if r.is_duplicate:
                        continue
                if not r.mapping_quality >= self.min_mapq:
                    continue
                if not (r.has_tag('BD') or r.has_tag('AD')):
                    continue
                if not clusters:
                    cluster = Cluster(shm_dir=self.shm_dir, max_proper_size=self.max_proper_pair_size)
                    cluster.append(r)
                    clusters.append(cluster)
                    continue
                if clusters[-1].read_is_compatible(r):
                    clusters[-1].append(r)
                elif len(clusters) >= 2 and r.reference_start == clusters[-1][-1].reference_start == clusters[-2][-1].reference_start:
                        clusters[-1].abnormal = True
                        clusters[-2].abnormal = True
                        # Could be X:22,432,984-22,433,240, a huge accumulation of fragments with rover homology.
                else:
                    cluster = Cluster(shm_dir=self.shm_dir, max_proper_size=self.max_proper_pair_size)
                    cluster.append(r)
                    clusters.append(cluster)
        logging.info('Found %d cluster on first pass (%s)', len(clusters), self.region or 0)
        if clusters:
            minimum_start = clusters[0].min
            maximum_end = clusters[-1].max
            cluster_density = len(clusters) / float(maximum_end - minimum_start)
            if cluster_density > 0.1:
                # every 10th nt a cluster, that should only happen on decoys.
                self.is_decoy = True
                logging.info('Skipping region with abnormally high cluster density (%s), probably a decoy (%s)',
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
                logging.info("Joining clusters (currently %d), round %i (%s)", cluster_length, i, self.region or 0)
                cluster_length = new_clusterlength
                for cluster in self.clusters:
                    cluster.join_adjacent(all_clusters=self.clusters)
                new_clusterlength = len(self.clusters)
        logging.info("Found %d cluster after first pass of cluster joining (%s).", new_clusterlength, self.region or 0)
        logging.info("Splitting cluster at polarity switches")
        for index, cluster in enumerate(self.clusters):
            new_clusters = cluster.split_cluster_at_polarity_switch()
            self._add_new_clusters(new_clusters, index)
        logging.info("After splitting at polarity switches we have %d cluster (was: %d) (%s)",
                     len(self.clusters),
                     new_clusterlength,
                     self.region or 0)
        logging.info("Checking cluster consistency (%s)", self.region or 0)
        for index, cluster in enumerate(self.clusters):
            new_clusters = cluster.check_cluster_consistency()
            self._add_new_clusters(new_clusters, index)
        logging.info("After splitting inconsistent clusters we have %d cluster", len(self.clusters))
        logging.info("Last pass of joining cluster (%s)", self.region or 0)
        for cluster in self.clusters:
            cluster.refine_members(self.assembly_realigner)
            cluster.join_adjacent(all_clusters=self.clusters)
        # We are done, we can give the clusters a numeric index, so that we can distribute the processing and recover the results
        logging.info("Found %d cluster overall (%s)", len(self.clusters), self.region or 0)
        [c.set_id(idx) for idx, c in enumerate(self.clusters)]

    def _add_new_clusters(self, new_clusters, index):
        current_index = index
        for cluster in new_clusters:
            if cluster:
                if current_index == index:
                    self.clusters[current_index] = cluster
                else:
                    self.clusters.insert(current_index, cluster)
                current_index += 1

    def collect_non_evidence(self):
        """Count reads that overlap cluster site but do not provide evidence for an insertion."""
        chunks = Chunks(clusters=self.clusters, header=self.header, input_path=self.input_path)
        logging.info("Collecting evidence (%s)", self.region or 0)
        for chunk in chunks.chunks:
            result = non_evidence(chunk)
            for index, evidence_against in result['against'].items():
                self.clusters[index].evidence_against = {r for l in evidence_against.values() for r in l}
                self.clusters[index].nref = len(evidence_against)
            for index, evidence_for in result['for'].items():
                for r, evidence in evidence_for.values():
                    if evidence == 'five_p':
                        self.clusters[index].evidence_for_five_p.add(r)
                    else:
                        self.clusters[index].evidence_for_three_p.add(r)

    def _create_contigs(self):
        futures = []
        for cluster in self.clusters:
            futures.append(self.tp.submit(cluster._make_contigs))
        wait(futures)

    def to_fasta(self):
        """Write supporting sequences to fasta file for detailed analysis."""
        logging.info("Writing contig fasta (%s)", self.region or 0)
        if self.output_fasta:
            self._create_contigs()
            with open(self.output_fasta, 'w') as out:
                for cluster in self.clusters:
                    for seq in cluster.to_fasta():
                        out.write(seq)

    def align_bwa(self):
        """Align cluster contigs or invidiual reads to a reference and write result into cluster."""
        logging.info("Aligning reads with BWA to describe cluster (%s)", self.region or 0)
        if self.output_fasta and (self.transposon_reference_fasta or self.transposon_bwa_index):
            bwa = Bwa(input_path=self.output_fasta,
                      bwa_index=self.transposon_bwa_index,
                      reference_fasta=self.transposon_reference_fasta,
                      threads=self.threads)
            for i, cluster in enumerate(self.clusters):
                description = bwa.desciption.get(i)
                if description:
                    cluster.insert_reference_name = description.pop(-1)
                cluster.feature_args = description

    def to_bam(self):
        """Write clusters of reads and include cluster number in CD tag."""
        logging.info("Writing clusters of reads (%s)", self.region or 0)
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


class Chunks(object):
    """Hold chunks of clusters."""

    def __init__(self, clusters, header, input_path):
        """Establish chunks of clusters for all clusters in `clusters`."""
        super(Chunks, self).__init__()
        self.clusters = clusters  # All clusters
        self.header = header
        self.input_path = input_path
        self.chunks = self._chunks()

    def _chunks(self):
        regions_d = {}
        for cluster in self.clusters:
            if cluster.abnormal:
                continue
            if cluster.tid not in regions_d:
                regions_d[cluster.tid] = []
            regions_d[cluster.tid].append(cluster)
        self.regions_d = regions_d  # Chromosome TID as key, clusters as values
        # stream over chromosomes and collect non support evidence
        allchunks = self._get_chunks()
        tasks = []
        for tid, chromosome_chunks in six.iteritems(allchunks):
            chromosome = self.header['SQ'][tid]['SN']
            input_path = self.input_path
            for chunk in chromosome_chunks:
                tasks.append({'chromosome': chromosome, 'input_path': input_path, 'chunk': chunk})
        return tasks

    def _get_chunks(self):
        """Return a list of chunks per chromosome."""
        chromosomes = {}
        for tid, cluster_list in six.iteritems(self.regions_d):
            chunks = self.chunk(cluster_list)
            chromosomes[tid] = ([chunk for chunk in chunks])
        return chromosomes

    @staticmethod
    def chunk(l, n=500):
        """Yield successive n-sized chunks from l."""
        for i in range(0, len(l), n):
            yield [chunk.serialize() for chunk in l[i:i + n]]
