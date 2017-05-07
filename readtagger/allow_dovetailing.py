import argparse
import logging

import pysam

from .bam_io import BamAlignmentWriter as Writer


logger = logging.getLogger(__name__)
logging.basicConfig(format='%(asctime)s %(name)s %(levelname)s - %(message)s', level=logging.DEBUG)


def get_max_proper_pair_size(alignment_file, reads_to_check=1000):
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
    for r in alignment_file:
        if r.is_proper_pair and not r.is_secondary and not r.is_supplementary:
            isize.append(abs(r.isize))
        if len(isize) == reads_to_check:
            alignment_file.reset()
            logger.info(msg, max(isize))
            return max(isize)
    alignment_file.reset()
    if isize:
        logger.info(msg, max(isize))
        return max(isize)
    else:
        logger.warn("Could not determine maximum allowed insert size for a proper pair. Are there any proper pairs in the input file?")
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
        logger.warn("Using default maximum insert size of %d" % default_max_proper_size)
        max_proper_size = default_max_proper_size
    if not read.is_proper_pair and not read.is_reverse == read.mate_is_reverse and read.reference_id == read.mrnm and abs(read.isize) <= max_proper_size:
        read.is_proper_pair = True
    return read


def main(args=None):
    """
    Main entrypoint for script.

    :param args:
    :return:
    """
    if not args:
        args = parse_args()
    with pysam.AlignmentFile(args.input_path) as input, Writer(args.output_path, template=input) as output:
        max_isize = get_max_proper_pair_size(input)
        [output.write(allow_dovetailing(read, max_proper_size=max_isize)) for read in input]


def parse_args():
    """
    Parse commandline arguments.

    :return: args
    :rtype argparse.ArgumentParser
    """
    from . import VERSION
    parser = argparse.ArgumentParser(description="Allow dovetailing.")
    parser.add_argument('-i', '--input_path', help="Input alignment file to manipulate", required=True)
    parser.add_argument('-o', '--output_path', help="Output alignment file", required=True)
    parser.add_argument('--version', action='version', version=VERSION)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
