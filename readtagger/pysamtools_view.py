from .bam_io import (
    BamAlignmentReader,
    BamAlignmentWriter
)


def view(input_bam, output_bam, region):
    """
    Extract reads from `input_bam` to `output_bam` for `region`.

    Adjusts coordinates so that the region start is 0 in the output bam file.
    """
    chromosome, start_end = region.split(':')
    start, end = start_end.split('-')
    start, end = int(start), int(end)
    with BamAlignmentReader(path=input_bam, sort_order='coordinate', region=region) as input:
        new_header = input.header.copy()
        new_header['SQ'] = [{'SN': chromosome, 'LN': end - start}]
        with BamAlignmentWriter(path=output_bam, header=new_header) as out:
            for read in input:
                read.next_reference_id = -1 if (read.next_reference_id != read.reference_id or read.next_reference_id == -1) else 0
                read.next_reference_start -= start
                if not start <= read.next_reference_start <= end:
                    read.mate_is_unmapped = True
                if read.next_reference_id < 0 or read.next_reference_start < 0 or read.mate_is_unmapped or read.is_unmapped:
                    read.mate_is_unmapped = True
                    read.next_reference_id = -1
                    read.next_reference_start = -1
                    read.is_proper_pair = False
                read.reference_id = 0
                read.reference_start -= start
                out.write(read)
