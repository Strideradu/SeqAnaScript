import argparse
import sys
from Bio import SeqIO

def rm_duplicate(args):
    records = SeqIO.parse(args.input, args.type)
    records = list(set(list(records)))
    SeqIO.write(records, args.output, args.type)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of input file", type=str)
    parser.add_argument("output", help="path of output file", type=str)
    parser.add_argument('--type', type=str, default='fasta', help='type of input "fasta" or "fastq"')

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    rm_duplicate(args)