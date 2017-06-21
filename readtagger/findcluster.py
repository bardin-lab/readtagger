import logging
import os
import shutil
import tempfile
from cached_property import cached_property
from concurrent.futures import (
    wait,
    ThreadPoolExecutor,
    ProcessPoolExecutor
)
import multiprocessing_logging
import pysam
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

logger = logging.getLogger(__name__)
logging.basicConfig(format='%(asctime)s %(name)s %(levelname)s - %(message)s', level=logging.DEBUG)
multiprocessing_logging.install_mp_handler()


class ClusterManager(object):
    """Coordinate multiple ClusterFinder objects when running in multiprocessing mode."""

    def __init__(self, **kwds):
        """Decide if passing kwds on to ClusterFinder or if splitting input file is required."""
        if kwds.get('max_proper_pair_size', 0) == 0:
            kwds['max_proper_pair_size'] = get_max_proper_pair_size(pysam.AlignmentFile(kwds['input_path']))
        if kwds['threads'] > 1:
            self.threads = kwds['threads']
            kwds['threads'] = 2
            self.kwds = kwds
            self.process_list = []
            self.process()
        else:
            # Delegate to clusterfinder
            ClusterFinder(**kwds)

    def process(self):
        """Process input bam in chunks."""
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            futures = []
            tempdir = tempfile.mkdtemp()
            chunks = split_locations_between_clusters(self.kwds['input_path'])
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
        self.merge_outputs()
        shutil.rmtree(tempdir, ignore_errors=True)

    def merge_outputs(self):
        """Merge outputs produced by working over smaller chunks with ClusterManager."""
        output_bam = self.kwds.get('output_bam')
        if output_bam:
            bam_files = [kwd['output_bam'] for kwd in self.process_list if kwd['output_bam']]
            merge_bam(bam_collection=bam_files, template_bam=bam_files[0], output_path=self.kwds['output_bam'])
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


def wrapper(kwds):
    """Launch ClusterFinder instances."""
    ClusterFinder(**kwds)


class ClusterFinder(object):
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
                 min_mapq=1,
                 max_clustersupport=800,
                 max_proper_pair_size=0,
                 remove_supplementary_without_primary=False,
                 region=None,
                 shm_dir=None):
        """
        Find readclusters in input_path file.

        This class assumes that all reads in input_path that should be clustered have either an 'AD' tag or a 'BD' tag.
        The initial limits of the cluster are defined by the maximum extent of overlapping reads, and each read that is added at the 3 prime end of
        the cluster will extend the cluster.
        The join_cluster method will then join clusters that overlap through their clipped sequences and cluster that can be assembled based on their proximity
        and the fact that they support the same same insertion (and can hence contribute to the same contig if assembled).
        """
        self.shm_dir = shm_dir
        self.region = region
        self.input_path = input_path
        self._sample_name = sample_name
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
        self._tempdir = tempfile.mkdtemp()
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
        self.cluster = self.find_cluster()
        self.clean_clusters()
        self.join_clusters()
        self.to_fasta()
        self.align_bwa()
        self.collect_non_evidence()
        self.to_bam()
        self.to_gff()
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def setup_bwa_indexes(self):
        """Handle setting up BWA indexes."""
        transposon_bwa_index = self.transposon_bwa_index
        genome_bwa_index = self.genome_bwa_index
        if not transposon_bwa_index and self.transposon_reference_fasta:
            transposon_bwa_index, _ = make_bwa_index(self.transposon_reference_fasta, dir=self._tempdir)
        if not genome_bwa_index and self.genome_reference_fasta:
            genome_bwa_index, _ = make_bwa_index(self.genome_reference_fasta, dir=self._tempdir)
        return transposon_bwa_index, genome_bwa_index

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

    def _remove_supplementary_without_primary(self):
        """Remove supplementary reads without primary alignments."""
        output_path = os.path.join(self._tempdir, 'clean.bam')
        discard_supplementary(input_path=self.input_path, output_path=output_path)
        self.input_path = output_path

    def find_cluster(self):
        """Find clusters by iterating over input_path and creating clusters if reads are disjointed."""
        if self.remove_supplementary_without_primary:
            self._remove_supplementary_without_primary()
        clusters = []
        with Reader(self.input_path, region=self.region, index=True) as reader:
            self.header = reader.header
            for r in reader:
                if not self.include_duplicates:
                    if r.is_duplicate:
                        continue
                if not r.mapping_quality > self.min_mapq:
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
                else:
                    cluster = Cluster(shm_dir=self.shm_dir, max_proper_size=self.max_proper_pair_size)
                    cluster.append(r)
                    clusters.append(cluster)
        logging.info('Found %d cluster', len(clusters))
        return clusters

    def clean_clusters(self):
        """Remove clusters that have more reads supporting an insertion than specified in self.max_clustersupport."""
        self.cluster = [c for c in self.cluster if not len(c.read_index) > self.max_clustersupport]

    def join_clusters(self):
        """Iterate over self.cluster and attempt to join consecutive clusters."""
        if len(self.cluster) > 1:
            cluster_length = len(self.cluster)
            new_clusterlength = 0
            while new_clusterlength != cluster_length:
                # Iterate until the length of the cluster doesn't change anymore.
                self._cache_join()
                cluster_length = new_clusterlength
                prev_cluster = self.cluster[0]
                for cluster in self.cluster[1:]:
                    if prev_cluster.can_join(cluster):
                        prev_cluster.extend(cluster)
                        self.cluster.remove(cluster)
                    else:
                        prev_cluster = cluster
                new_clusterlength = len(self.cluster)
        for index, cluster in enumerate(self.cluster):
            cluster_a, cluster_b = cluster.split_cluster_at_polarity_switch()
            if cluster_a:
                self.cluster[index] = cluster_a
            if cluster_b:
                self.cluster.insert(index + 1, cluster_b)
        for cluster in self.cluster:
            cluster.refine_members(self.assembly_realigner)
        # We are done, we can give the clusters a numeric index, so that we can distribute the processing and recover the results
        [c.set_id(i) for i, c in enumerate(self.cluster)]

    def _cache_join(self):
        futures = []
        prev_cluster = self.cluster[0]
        for cluster in self.cluster[1:]:
            futures.append(self.tp.submit(prev_cluster.can_join, cluster))
            prev_cluster = cluster
        wait(futures)

    def collect_non_evidence(self):
        """Count reads that overlap cluster site but do not provide evidence for an insertion."""
        chunks = Chunks(clusters=self.cluster, header=self.header, input_path=self.input_path)
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            r = executor.map(non_evidence, chunks.chunks)
            for result in r:
                for index, nref in result.items():
                    self.cluster[index].nref = nref

    def _create_contigs(self):
        futures = []
        for cluster in self.cluster:
            futures.append(self.tp.submit(cluster._make_contigs))
        wait(futures)

    def to_fasta(self):
        """Write supporting sequences to fasta file for detailed analysis."""
        if self.output_fasta:
            self._create_contigs()
            with open(self.output_fasta, 'w') as out:
                for cluster in self.cluster:
                    for seq in cluster.to_fasta():
                        out.write(seq)

    def align_bwa(self):
        """Align cluster contigs or invidiual reads to a reference and write result into cluster."""
        if self.output_fasta and self.transposon_reference_fasta:
            bwa = Bwa(input_path=self.output_fasta,
                      bwa_index=self.transposon_bwa_index,
                      reference_fasta=self.transposon_reference_fasta,
                      threads=self.threads)
            for i, cluster in enumerate(self.cluster):
                description = bwa.desciption.get(i)
                if description:
                    cluster.reference_name = description.pop(-1)
                cluster.feature_args = description

    def to_bam(self):
        """Write clusters of reads and include cluster number in CD tag."""
        if self.output_bam:
            with Writer(self.output_bam, header=self.header) as writer:
                for i, cluster in enumerate(self.cluster):
                    for r in cluster:
                        r.set_tag('CD', i)
                        writer.write(r)
            # Because the cluster splitting doesn't necessarily conserve order we need to sort again.
            sort_bam(inpath=self.output_bam, output=self.output_bam, sort_order='coordinate', threads=self.threads)

    def to_gff(self):
        """Write clusters as GFF file."""
        if self.output_gff:
            write_cluster(clusters=self.cluster,
                          header=self.header,
                          output_path=self.output_gff,
                          sample=self.sample_name,
                          threads=self.threads)
            sort_gff(self.output_gff, output_path=self.output_gff)


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
