import os
from collections import OrderedDict
from functools import partial

from concurrent.futures import ThreadPoolExecutor
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from .utils import sort


def write_gff_cluster(clusters, header, output_path, sample_name='sample', threads=1):
    """Write clusters as GFF entries."""
    with open(output_path, "w") as out_handle:
        tp = ThreadPoolExecutor(threads)
        futures = []
        records = OrderedDict((tid, SeqRecord(Seq(""), sn)) for tid, sn in enumerate(header.references))
        for i, cluster in enumerate(clusters):
            if not cluster.exclude:
                func = partial(get_feature, cluster, sample_name, i)
                futures.append(tp.submit(func))
        for future in futures:
            tid, feature = future.result()
            records[tid].features.append(feature)
        GFF.write(records.values(), out_handle)
        tp.shutdown(wait=True)


def merge_gff_files(gff_files, output_path, sort=True):
    """Merge multiple GFF files."""
    wrote_header = False
    with open(output_path, 'w') as gff_writer:
        for output_piece in gff_files:
            if os.path.exists(output_piece):
                with open(output_piece) as piece:
                    for line in piece:
                        if line.startswith('#') and wrote_header:
                            continue
                        gff_writer.write(line)
            wrote_header = True
    if sort:
        sort_gff(input_path=output_path, output_path=output_path)


def get_feature(cluster, sample, i):
    """Turn a cluster into a biopython SeqFeature."""
    qualifiers = OrderedDict(ID="%s_%d" % (sample, i))
    for attr in cluster.exportable:
        qualifiers[attr] = getattr(cluster, attr.lower())
    feature = SeqFeature(FeatureLocation(cluster.start, cluster.end), type=cluster.type, strand=1, qualifiers=qualifiers)
    if cluster.feature_args:
        seqfeature_args = cluster.feature_args
        base_args = {'location': FeatureLocation(cluster.start, cluster.end), 'strand': 1}
        subfeatures = []
        for feature_args in seqfeature_args:
            feature_args = feature_args.to_feature_args()
            args = base_args.copy()
            args.update(feature_args)
            subfeatures.append(SeqFeature(**args))
        feature.sub_features = subfeatures
    return cluster.tid, feature


def sort_gff(input_path, output_path):
    """Sort GFF file."""
    return sort(input_path=input_path, output_path=output_path, sort_cmd="sort -k 1,1 -k4,4n")
