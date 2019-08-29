import logging

import pysam

from .bam_io import BamAlignmentWriter as Writer

logger = logging.getLogger(__name__)


def get_max_proper_pair_size(path, reads_to_check=1000, region=None):
    """
    Iterate over the first 1000 properly paired records in alignment_file and get the maximum valid isize for a proper pair.

    :param alignment_file: pysam.AlignmentFile
    :type alignment_file: pysam.AlignmentFile
    :param reads_to_check: number of reads to check
    :type int
    :rtype int
    """
    isize = []
    msg = "Maximum insert size for a proper pair is %s"
    with pysam.AlignmentFile(path) as alignment_file:
        if region:
            iter_reads = alignment_file.fetch(region=region)
        else:
            iter_reads = alignment_file
        for r in iter_reads:
            if r.is_proper_pair and not r.is_secondary and not r.is_supplementary:
                isize.append(abs(r.isize))
            if len(isize) == reads_to_check:
                alignment_file.reset()
                logger.info(msg, max(isize))
                return max(isize)
    if isize:
        logger.info(msg, max(isize))
        return max(isize)
    else:
        logger.warning("Could not determine maximum allowed insert size for a proper pair. Are there any proper pairs in the input file?")
        return None


def allow_dovetailing(read, max_proper_size=351, default_max_proper_size=351):
    """
    Manipulate is_proper_pair tag to allow dovetailing of reads.

    Precondition is read and mate have the same reference id, are within the maximum proper pair distance
    and are either in FR or RF orientation.
    :param read: aligned segment of pysam.AlignmentFile
    :type read: pysam.AlignedSegment
    :rtype pysam.AlignedSegment
    """
    if max_proper_size is None:
        logger.warning("Using default maximum insert size of %d" % default_max_proper_size)
        max_proper_size = default_max_proper_size
    if not read.is_proper_pair and not read.is_reverse == read.mate_is_reverse and read.reference_id == read.mrnm and abs(read.isize) <= max_proper_size:
        read.is_proper_pair = True
    return read


def process(input_path, output_path):
    """
    Run script.

    :param args:
    :return:
    """
    max_isize = get_max_proper_pair_size(input_path)
    with pysam.AlignmentFile(input_path) as input, Writer(output_path, template=input) as output:
        [output.write(allow_dovetailing(read, max_proper_size=max_isize)) for read in input]
