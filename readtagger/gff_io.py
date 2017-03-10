import collections
from functools import partial
import yaml
from concurrent.futures import ThreadPoolExecutor
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from .hits import BlastProcessor


def write_cluster(clusters, header, output_path, reference_fasta=None, blastdb=None, sample='sample', threads=1):
    """Write clusters as GFF entries."""
    with open(output_path, "w") as out_handle:
        if blastdb or reference_fasta:
            blast = BlastProcessor(reference_fasta=reference_fasta, blastdb=blastdb)
        else:
            blast = None
        tp = ThreadPoolExecutor(threads)
        futures = []
        records = {tid: SeqRecord(Seq(""), header['SQ'][tid]['SN']) for tid in range(len(header['SQ']))}
        for i, cluster in enumerate(clusters):
            func = partial(get_feature, cluster, sample, i, blast)
            futures.append(tp.submit(func))
        for future in futures:
            tid, feature = future.result()
            records[tid].features.append(feature)
        GFF.write(records.values(), out_handle)
        tp.shutdown(wait=True)


def get_feature(cluster, sample, i, blast=None):
    """Turn a cluster into a biopython SeqFeature."""
    qualifiers = {"source": "findcluster",
                  "score": cluster.score,
                  "left_support": cluster.left_support,
                  "right_support": cluster.right_support,
                  "left_insert": cluster.left_contigs,
                  "right_insert": cluster.right_contigs,
                  "ID": "%s_%d" % (sample, i),
                  "valid_TSD": cluster.valid_tsd}
    if blast:
        left_description, right_description, common_inserts = blast.blast_cluster(cluster)
        qualifiers['left_description'] = yaml.dump(convert(left_description), default_flow_style=False).replace('\n', ' ').replace('-', '')
        qualifiers['right_description'] = yaml.dump(convert(right_description), default_flow_style=False).replace('\n', ' ').replace('-', '')
        qualifiers['common_inserts'] = yaml.dump(list(convert(common_inserts)), default_flow_style=False).replace('\n', ' ').replace('-', '')

    return cluster.tid, SeqFeature(FeatureLocation(cluster.start, cluster.end), type="TE", strand=1, qualifiers=qualifiers)

def convert(data):
    """Convert unicode to bytestring"""
    if isinstance(data, basestring):
        return str(data)
    elif isinstance(data, collections.Mapping):
        return dict(map(convert, data.iteritems()))
    elif isinstance(data, collections.Iterable):
        return type(data)(map(convert, data))
    else:
        return data
