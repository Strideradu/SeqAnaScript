import argparse
import sys
import os
from Bio import SeqIO
from collections import defaultdict
import SAMparser
import numpy as np
import dill as pickle

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of input file", type=str)
    # parser.add_argument("fasta", help="path of fasta file", type=str)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    aligns = defaultdict(lambda: defaultdict(list))
    with open(args.input) as f:
        for line in f:
            if line[0] != "@":
                record = SAMparser.text_to_sam(line)

                id = record.rname
                tlen = record.tlen

                if id != '*':
                    aligns[record.qname][tlen].append(record)

    pos_aligns = []
    for key, tlen_dict in aligns.items():
        pos = []
        texts = []
        is_chimera = False
        if len(tlen_dict) > 2:
            is_chimera = True

        for tlen, records in tlen_dict.items():
            if tlen > 0:
                if -tlen not in tlen_dict:
                    is_chimera = True

            for record in records:
                pos.append(record.pos)
                texts.append(record.text)
                if record.sa:
                    if record.sa_rname == record.rname:
                        pos.append(record.sa_pos)

        min_pos = np.min(pos)
        max_pos = np.max(pos)

        if max_pos - min_pos >= 1000:
            pos_aligns.append((min_pos, min_pos, max_pos))
            pos_aligns.append((max_pos, min_pos, max_pos))

    clusters = [[0, 0]]
    last = 0
    pos_aligns.sort()
    for aligns in pos_aligns:
        align_pos, _, _ = aligns
        if align_pos >= last + 150:
            clusters.append([align_pos, align_pos + 150])
            last = align_pos + 150

        else:
            clusters[-1][1] = max(clusters[-1][1], align_pos + 150)
            last = clusters[-1][1]

    print(clusters)
