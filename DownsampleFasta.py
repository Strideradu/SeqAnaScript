import argparse
import sys
from Bio import SeqIO
import numpy as np
from scipy.stats.kde import gaussian_kde
import random


def get_pdf(lens):
    kde_pdf = gaussian_kde(lens)
    return kde_pdf


def get_coverage(lens, genome_len):
    return np.sum(lens) / genome_len


def random_downsample(path, output, target, genome_len):
    length = []
    with open(path) as f:
        records = list(SeqIO.parse(f, 'fasta'))
        for record in records:
            length.append(len(record.seq))

    cov = get_coverage(length, genome_len)
    pdf = get_pdf(length)

    ratio = target/cov
    p_0 = 1/(np.max(length) - np.min(length))
    downsample = []
    down_len = []
    for record in records:
        read_len = len(record.seq)
        rd = pdf.pdf(read_len)/p_0 * ratio
        if random.random() <= rd:
            downsample.append(record)
            down_len.append(read_len)

    print("After downsample, the total read length is {}".format(np.sum(down_len)))
    print("Coverage is {}".format(np.sum(down_len)/genome_len))

    SeqIO.write(downsample, output, 'fasta')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of input file", type=str)
    parser.add_argument("output", help="path of output file", type=str)
    parser.add_argument("target", help="target coverage", type=int)
    parser.add_argument("--genome-path", help="path of genome file", type=str, default=None)
    parser.add_argument("--genome-len", help="length of genome", type=int, default=None)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    if not(args.genome_path) and not(args.genome_len):
        print("Must provide genome file path or genome length")
        parser.print_help()
        sys.exit(1)

    if args.genome_path:
        record = SeqIO.read(args.genome_path, 'fasta')
        genome_len = len(record.seq)
    else:
        genome_len = args.genome_len

    random_downsample(args.input, args.output,args.target, genome_len)