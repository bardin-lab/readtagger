import argparse
import pysam
from .bam_io import BamAlignmentWriter as Writer
from .readtagger import SamAnnotator
from .readtagger import __VERSION__


def main(args=None):
    """
    Main entrypoint for script.

    :param args:
    :return:
    """
    if not args:
        args = parse_args()
    with pysam.AlignmentFile(args.input_path) as input, Writer(args.output_path, template=input) as output:
        max_isize = SamAnnotator.get_max_proper_pair_size(input)
        [output.write(SamAnnotator.allow_dovetailing(read, max_proper_size=max_isize)) for read in input]


def parse_args():
    """
    Parse commandline arguments.

    :return: args
    :rtype argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(description="Allow dovetailing.")
    parser.add_argument('-i', '--input_path', help="Input alignment file to manipulate", required=True)
    parser.add_argument('-o', '--output_path', help="Output alignment file", required=True)
    parser.add_argument('--version', action='version', version=__VERSION__)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
