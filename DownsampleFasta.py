import argparse
import random
import sys

import numpy as np
from Bio import SeqIO
from scipy.stats.kde import gaussian_kde


def weighted_sample(population, weights, k):
    """
    Alternative way to previous implementation.

    This function draws a random sample of length k
    from the sequence 'population' according to the
    list of weights
    """
    sample = set()
    population = list(population)
    weights = list(weights)
    while len(sample) < k:
        choices = random.choices(population, weights)
        for choice in choices:
            if choice not in sample:
                sample.add(choice)
    return list(sample)


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
            seq_len = len(record.seq)
            length.append(seq_len)

    if genome_len == 1:
        ratio = target
    else:
        cov = get_coverage(length, genome_len)
        ratio = target / cov
    pdf = get_pdf(length)

    size = int(len(length) * ratio)
    weight = [pdf.pdf(x) for x in length]
    downsample = []
    downsample = weighted_sample(records, weight, size)
    down_len = [len(x.seq) for x in downsample]

    print("After downsample, the total read length is {}".format(np.sum(down_len)))
    print("Coverage is {}".format(np.sum(down_len) / genome_len))

    SeqIO.write(downsample, output, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of input file", type=str)
    parser.add_argument("output", help="path of output file", type=str)
    parser.add_argument("target", help="target coverage", type=int)
    parser.add_argument("--genome-path", help="path of genome file", type=str, default=None)
    parser.add_argument("--genome-len", help="length of genome", type=int, default=None)
    parser.add_argument("--ratio", help="ratio to sample", type=float, default=None)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    if not (args.genome_path) and not (args.genome_len) and not (args.ratio):
        print("Must provide genome file path or genome length or sample ratio")
        parser.print_help()
        sys.exit(1)

    if args.ratio:
        genome_len = 1
        args.target = args.ratio
    elif args.genome_path:
        record = SeqIO.read(args.genome_path, 'fasta')
        genome_len = len(record.seq)
    else:
        genome_len = args.genome_len

    random_downsample(args.input, args.output, args.target, genome_len)
