import os
from collections import (
    namedtuple,
    OrderedDict
)
import pysam
from edlib import align

gff_record = namedtuple('GFF_record', 'seqid source type start end score strand phase attributes')
comparison = namedtuple('Comparison', 'putative_event treatment_events control_events')

SCAN_SOFT_CLIP_REGION = 15


def get_tabix_file(path, preset='gff'):
    """Return a TabixFile instance for a file at `path`."""
    if not os.path.exists("%s.gz.tbi" % path):
        pysam.tabix_index(path, preset=preset, keep_original=True)
    tabix_file = pysam.TabixFile("%s.gz" % path)
    return tabix_file


def to_gff_attributes(s):
    """Convert a GFF attribute string to an OrderedDict."""
    attr_d = OrderedDict()
    attributes = s.split(';')
    for attribute in attributes:
        k, v = attribute.split('=')
        if ',' in v:
            v = v.split(',')
        attr_d[k] = v
    return attr_d


def to_gff_record(r):
    """Convert a string representation of a GFF record to a named tuple."""
    fields = r.split('\t')
    fields[3], fields[4] = int(fields[3]), int(fields[4])
    fields[8] = to_gff_attributes(fields[8])
    return gff_record(*fields)


def gff_record_to_string(gff_record):
    """Convert a `GFF_record` namedtuple to a string representation."""
    fields = list(gff_record)
    fields[3], fields[4] = str(fields[3]), str(fields[4])
    attributes = []
    for k, v in fields[8].items():
        if isinstance(v, list):
            v = ",".join(v)
        attributes.append("%s=%s" % (k, v))
    fields[8] = ";".join(attributes)
    return "%s\n" % "\t".join(fields)


def write_gff_records(records, path):
    """Write GFF_record in records iterable to `path`."""
    with open(path, 'w') as out:
        for r in records:
            out.write(gff_record_to_string(r))


def fetch_records(tabixfile, region):
    """Wrap TabixFile.fetch in try/except clause to ignore region errors."""
    try:
        return tabixfile.fetch(reference=region.seqid,
                               start=int(region.start) - SCAN_SOFT_CLIP_REGION,
                               end=int(region.end) + SCAN_SOFT_CLIP_REGION)

    except ValueError:
        # Happens if region is not present in tabix file, that's fine and expected
        return []


def fill_comparison(putative, controls, treatment):
    """Convert tabix record to a Comparison named tuple."""
    for putative_event in putative.fetch():
        putative_record = to_gff_record(putative_event)
        putative_record_complements = [to_gff_record(r) for r in fetch_records(treatment, putative_record)]
        control_records = []
        for control in controls:
            for r in fetch_records(control, putative_record):
                control_records.append(to_gff_record(r))
        yield (comparison(putative_record, putative_record_complements, control_records))


def annotate_with_overlapping_insertions(putative_record, control_records):
    """If overlapping insertion of the same type is found list them in the overlaps attribute."""
    putative_insert_reference = putative_record.attributes.get('insert_reference_name')
    overlapping_ids = []
    for control_record in control_records:
        if control_record.source == 'findcluster' and control_record.attributes.get('insert_reference_name') == putative_insert_reference:
            overlapping_ids.append(control_record.attributes.get('ID'))
    putative_record.attributes['overlaps'] = ",".join(overlapping_ids)
    return putative_record


def filter_putative_insertions(putative, treatment, controls, min_length, output_discarded_records=True):
    """Remove insertions that are based on clipped sequences which are also present in the control."""
    for (putative_record, putative_complements, control_records) in fill_comparison(putative, controls, treatment):
        valid_record = True
        clips = []
        control_clips = []
        putative_record = annotate_with_overlapping_insertions(putative_record, control_records)
        softclip_clusters = putative_record.attributes.get('softclip_clusters')
        if not softclip_clusters:
            yield putative_record
            continue
        for putative_complement in putative_complements:
            if putative_complement.type in ('3p_clip', '5p_clip') and putative_complement.attributes['ID'] in softclip_clusters:
                # We only verify clipped sequences that are part of the insertion to verify
                clips.append(putative_complement)
        for control_record in control_records:
            if control_record.type in ('3p_clip', '5p_clip'):
                control_clips.append(control_record)
        for c in control_clips:
            for t in clips:
                if c.type == t.type:
                    if abs(c.start - t.start) < SCAN_SOFT_CLIP_REGION or abs(c.end - t.end) < SCAN_SOFT_CLIP_REGION:
                        c_seq = c.attributes.get('consensus')
                        t_seq = t.attributes.get('consensus')
                        if sequences_match(seq1=c_seq, seq2=t_seq, compare=c.type, min_length=min_length):
                            valid_record = False
                            break
            if valid_record is False:
                break
        if valid_record:
            yield putative_record
        elif output_discarded_records:
            putative_record.attributes['FAIL'] = 'clip_seq_matches_%s' % c.attributes.get('ID')
            yield putative_record


def sequences_match(seq1, seq2, compare='5p_clip', min_length=4):
    """
    Compare 2 sequences and determine whether they are likely the same.

    >>> s1, s2 = 'AAAAAAAAAA', 'AAAAAAAAATGC'
    >>> sequences_match(s1, s2, compare='3p_clip')
    True
    >>> sequences_match('TTGAAAAAAA', s2, compare='3p_clip')
    False
    >>> sequences_match(s1[::-1], s2[::-1], compare='5p_clip')
    True
    >>> sequences_match('A', 'A')
    False
    """
    if compare == '3p_clip':
        seq1 = seq1[:10]
        seq2 = seq2[:10]
    else:
        seq1 = seq1[-10:]
        seq2 = seq2[-10:]
    if len(seq1) > min_length and len(seq2) > min_length:
        if len(seq1) >= 10 and len(seq2) >= 10:
            if align(seq1, seq2)['editDistance'] < 2:
                return True
        elif compare == '3p_clip':
            return seq1.startswith(seq2) or seq2.startswith(seq1)
        else:
            return seq1.endswith(seq2) or seq2.endswith(seq1)
    return False


def confirm_insertions(putative_insertions_path, all_treatments_path, all_controls_paths, output_path, min_length=4, output_discarded_records=True):
    """Confirm putative insertions by absence of softclipping patterns in a somatic control file."""
    putative = get_tabix_file(putative_insertions_path)
    treatment = get_tabix_file(all_treatments_path)
    controls = [get_tabix_file(p) for p in all_controls_paths if not p == all_treatments_path]
    gff_record_iterator = filter_putative_insertions(putative=putative,
                                                     treatment=treatment,
                                                     controls=controls,
                                                     min_length=min_length,
                                                     output_discarded_records=output_discarded_records)
    write_gff_records(gff_record_iterator, path=output_path)
