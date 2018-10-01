import argparse
import sys
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def to_original_seq(args):
    seqs = []
    for multiple_alignment in AlignIO.parse(args.input, "maf"):
        multiple_alignment = list(multiple_alignment)
        corrected_str = str(multiple_alignment[0].seq).replace("-", "")
        corrected_seq = Seq(corrected_str)
        record = SeqRecord(corrected_seq, id=multiple_alignment[1].id, name="", description="")
        seqs.append(record)

    with open(args.output, "w") as output_handle:
        SeqIO.write(seqs, output_handle, "fasta")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of input file", type=str)
    parser.add_argument("output", help="path of output file", type=str)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    to_original_seq(args)
