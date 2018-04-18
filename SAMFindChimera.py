import argparse
import sys
import os
from Bio import SeqIO
from collections import defaultdict
import SAMparser
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of input file", type=str)
    parser.add_argument("fasta", help="path of fasta file", type=str)
    parser.add_argument("output", help="path of output file", type=str)
    parser.add_argument("id_list", help="path of id list file", type=str)

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

                if id!='*':
                    aligns[id][tlen].append(record)

    records = SeqIO.parse(args.fasta, 'fasta')
    for record in records:
        if record.id == '277694':
            table_size = len(record.seq)//500 + 1

    count = np.zeros((table_size, table_size))
    sin_count = np.zeros((table_size, table_size))

    with open(args.output,'w') as f1:
        with open(args.id_list, 'w') as f2:
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
                            pos.append(record.sa_pos)
                            min_pos = min(record.pos, record.sa_pos)
                            max_pos = max(record.pos, record.sa_pos)

                            if max_pos - min_pos >= 1000:
                                sin_count[min_pos// 500 + 1][max_pos// 500 + 1] += 1

                min_pos = np.min(pos)
                max_pos = np.max(pos)

                if max_pos - min_pos >= 1000:
                    count[min_pos//500 + 1][max_pos// 500 + 1] += 1

                if is_chimera:
                    print(key, file=f2)
                    for text in texts:
                        print(text, file=f1)

    imgplot = plt.imshow(count[:, :, 0])
    plt.savefig("matrix_count.png")

    plt.figure()
    imgplot = plt.imshow(sin_count[:, :, 0])
    plt.savefig("sin_matrix_count.png")


