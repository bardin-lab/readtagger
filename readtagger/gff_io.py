import os
import subprocess
import tempfile
from functools import partial

from concurrent.futures import ThreadPoolExecutor
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def write_cluster(clusters, header, output_path, sample='sample', threads=1):
    """Write clusters as GFF entries."""
    with open(output_path, "w") as out_handle:
        tp = ThreadPoolExecutor(threads)
        futures = []
        records = {tid: SeqRecord(Seq(""), header['SQ'][tid]['SN']) for tid in range(len(header['SQ']))}
        for i, cluster in enumerate(clusters):
            func = partial(get_feature, cluster, sample, i)
            futures.append(tp.submit(func))
        for future in futures:
            tid, feature = future.result()
            records[tid].features.append(feature)
        GFF.write(records.values(), out_handle)
        tp.shutdown(wait=True)


def get_feature(cluster, sample, i):
    """Turn a cluster into a biopython SeqFeature."""
    genotype = cluster.genotype_likelihood()
    qualifiers = {"source": "findcluster",
                  "score": cluster.score,
                  "left_support": cluster.left_support,
                  "left_mate_support": cluster.left_mate_count,
                  "right_support": cluster.right_support,
                  "right_mate_support": cluster.right_mate_count,
                  "non_support": genotype.nref,
                  "genotype": genotype.genotype,
                  "genotype_likelihoods": [genotype.reference, genotype.heterozygous, genotype.homozygous],
                  "left_insert": [v for pair in enumerate(cluster.left_contigs) for v in pair],
                  "right_insert": [v for pair in enumerate(cluster.right_contigs) for v in pair],
                  "ID": "%s_%d" % (sample, i),
                  "valid_TSD": cluster.valid_tsd}
    if cluster.reference_name:
        qualifiers['insert_reference'] = cluster.reference_name
    feature = SeqFeature(FeatureLocation(cluster.start, cluster.end), type=cluster.reference_name or "TE", strand=1, qualifiers=qualifiers)
    if cluster.feature_args:
        seqfeature_args = cluster.feature_args
        base_args = {'location': FeatureLocation(cluster.start, cluster.end), 'strand': 1}
        subfeatures = []
        for feature_args in seqfeature_args:
            args = base_args.copy()
            args.update(feature_args)
            subfeatures.append(SeqFeature(**args))
        feature.sub_features = subfeatures
    return cluster.tid, feature


def sort_gff(input_path, output_path):
    """Sort gff file at input path."""
    fd, tmp = tempfile.mkstemp()
    header_lines = []
    with open(input_path) as gff_in, open(tmp, 'w') as out:
        for line in gff_in:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                out.write(line)
        out.close()
    with open(output_path, 'w') as out:
        p = subprocess.Popen(['sort', '-k', '1,1', '-k4,4n', tmp], stdout=subprocess.PIPE, close_fds=True)
        out.write("\n".join(header_lines))
        for line in p.stdout:
            out.write(line.decode())
    os.close(fd)
