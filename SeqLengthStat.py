import argparse
import sys
from Bio import SeqIO
import numpy as np

def length_stat(path, type):
    length = []
    if type == 'fastq':
        errors = [] # expected wrong bases in each read
    with open(path) as f:
        records = SeqIO.parse(f, type)
        for record in records:
            length.append(len(record.seq))
            if type == 'fastq':
                q = np.array(record.letter_annotations["phred_quality"])
                pe = 10**(-0.1 * q)
                errors.append(length[-1] * pe)

    length = np.array(length)
    count = len(length)
    max_len = np.max(length)
    min_len = np.min(length)
    mean = np.mean(length)
    median = np.median(length)
    total = np.sum(length)
    print("There are {} reads".format(count))
    print("In total, there are {} bases".format(total))
    print("The mean of reads is {}".format(mean))
    print("The median of reads is {}".format(median))
    print("longest length is {} and the shortest length is {}".format(max_len, min_len))
    if type == 'fastq':
        print('Estimated error rates is {}'.format(np.sum(errors)/total))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of input file", type=str)
    parser.add_argument('--type', type=str, default='fasta', help='type of input "fasta" or "fastq"')

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    length_stat(args.input, args.type)


