from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def write_cluster(clusters, header, output_path, sample='sample'):
    """Write clusters as GFF entries."""
    with open(output_path, "w") as out_handle:
        for i, cluster in enumerate(clusters):
            record = get_record(header=header, cluster=cluster, sample=sample, i=i)
            GFF.write([record], out_handle)


def get_record(header, cluster, sample, i):
    """Turn a cluster into a biopython SeqFeature."""
    tid = cluster[0].tid
    chromosome = header['SQ'][tid]['SN']
    rec = SeqRecord(Seq(""), chromosome)  # seems a bit silly, but OK ...
    left_sequences = list(cluster.clustertag.left_sequences.keys())
    right_sequences = list(cluster.clustertag.right_sequences.keys())
    all_sequences = left_sequences + right_sequences
    if cluster.clustertag.left_insert:
        left_contigs = (contig.sequence for contig in cluster.clustertag.left_insert.assembly.contigs)
    else:
        left_contigs = ()
    left_contigs = [v for pair in enumerate(left_contigs) for v in pair]
    if cluster.clustertag.right_insert:
        right_contigs = (contig.sequence for contig in cluster.clustertag.right_insert.assembly.contigs)
    else:
        right_contigs = ()
    right_contigs = [v for pair in enumerate(right_contigs) for v in pair]
    qualifiers = {"source": "findcluster",
                  "score": len(all_sequences),
                  "left_support": len(left_sequences),
                  "right_support": len(right_sequences),
                  "left_insert": left_contigs,
                  "right_insert": right_contigs,
                  "ID": "%s_%d" % (sample, i),
                  "valid_TSD": cluster.clustertag.tsd.is_valid}
    start = cluster.clustertag.five_p_breakpoint
    end = cluster.clustertag.three_p_breakpoint
    if not start:
        start = end
    if not end:
        end = start
    if start > end:
        end, start = start, end
    if start == end:
        end += 1
    top_feature = SeqFeature(FeatureLocation(start, end), type="TE", strand=1, qualifiers=qualifiers)
    rec.features = [top_feature]
    return rec
