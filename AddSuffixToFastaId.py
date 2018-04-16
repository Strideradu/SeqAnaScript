import argparse
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("input", help="path of input file", type=str)
parser.add_argument("output", help="path of output file", type=str)
parser.add_argument("suffix", help="suffix", type=str)

try:
    args = parser.parse_args()

except:
    parser.print_help()
    sys.exit(1)

out_records = []
with open(args.input) as f:
    records = list(SeqIO.parse(f, 'fasta'))
    for record in records:
        record.id = record.id + "_" + args.suffix
        out_records.append(record)

SeqIO.write(out_records, args.output, 'fasta')