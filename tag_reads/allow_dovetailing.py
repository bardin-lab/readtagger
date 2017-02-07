import argparse
import pysam
from tag_reads import SamAnnotator


def main():
    args = parse_args()
    alignment_file = pysam.AlignmentFile(args.input_path)
    output = pysam.AlignmentFile(args.output_path, 'wb', template=alignment_file)
    max_isize = SamAnnotator.get_max_proper_pair_size(alignment_file)
    [output.write(SamAnnotator.allow_dovetailing(read, max_proper_size=max_isize)) for read in alignment_file]


def parse_args():
    parser = argparse.ArgumentParser(description="Allow dovetailing.")
    parser.add_argument('-i', '--input_path', help="Input alignment file to manipulate", required=True)
    parser.add_argument('-o', '--output_path', help="Output alignment file", required=True)
    parser.parse_args()

if __name__ == '__main__':
    main()