# Combine all fasta files under one folder to single multi-fasta file
import argparse
import sys
import os
from Bio import SeqIO


def combine_files(args):
    """
    combine all fasta files
    :param path:
    :return:
    """
    files = os.listdir(args.input)
    output = []
    for file in files:
        try:
            records = SeqIO.parse(file, format=args.format)
            output += list(records)
        except:
            continue

    if len(output) > 0:
        SeqIO.write(output, args.output, format=args.format)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of the folder", type=str)
    parser.add_argument("output", help="path of output file", type=str)
    parser.add_argument("--format", help="combine what format file", default='fasta', type=str)

    try:
        args = parser.parse_args()
        assert os.path.isdir(args.input) is True

    except:
        parser.print_help()
        sys.exit(1)
