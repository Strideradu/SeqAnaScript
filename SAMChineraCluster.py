import argparse
import sys
import os
from Bio import SeqIO
from collections import defaultdict
import SAMparser
import numpy as np
import dill as pickle
import intervaltree


def build_intervaltree(input):
    gene_pos = {}
    tree = intervaltree.IntervalTree()
    with open(input) as f:
        for line in f:
            if line[0] != "#":
                sp = line.strip().split()
                if sp[3] == "277694":
                    start = min(int(sp[5]), int(sp[4]))
                    end = max(int(sp[5]), int(sp[4]))
                    name = sp[8]
                    tree[start:end + 1] = name
                    gene_pos[name] = (start, end + 1)

    return tree, gene_pos


def check_annotation(tree, gene_pos, count, clusters, output):
    with open(output, "w") as fout:
        for ct in count:

            cluster1 = ct[1]
            cluster2 = ct[2]
            freq = ct[0]
            res1 = tree[clusters[cluster1][0]:clusters[cluster1][1]]
            res2 = tree[clusters[cluster2][0]:clusters[cluster2][1]]

            if (len(res1) == 0 and len(res2) == 0):
                print(
                    "f\t{}\talign1\t{}\t{}\t{}\t{}\t{}\talign2\t{}\t{}\t{}\t{}\t{}".format(freq, clusters[cluster1][0],
                                                                                           clusters[cluster1][1], "*",
                                                                                           "*", "*",
                                                                                           clusters[cluster2][0],
                                                                                           clusters[cluster2][1], "*",
                                                                                           "*", "*"), file=fout)
                continue

            elif (len(res1) != 0 and len(res2) != 0):
                for interval1 in res1:
                    for interval2 in res2:
                        name1 = interval1.data
                        name2 = interval2.data
                        print("f\t{}\talign1\t{}\t{}\t{}\t{}\t{}\talign2\t{}\t{}\t{}\t{}\t{}".format(freq,
                                                                                                     clusters[cluster1][
                                                                                                         0],
                                                                                                     clusters[cluster1][
                                                                                                         1], name1,
                                                                                                     gene_pos[name1][0],
                                                                                                     gene_pos[name1][1],
                                                                                                     clusters[cluster2][
                                                                                                         0],
                                                                                                     clusters[cluster2][
                                                                                                         1], name2,
                                                                                                     gene_pos[name2][0],
                                                                                                     gene_pos[name2][
                                                                                                         1]), file=fout)

            else:
                if len(res1) != 0:
                    for interval in res1:
                        name = interval.data

                        print("t\t{}\talign1\t{}\t{}\t{}\t{}\t{}\talign2\t{}\t{}\t{}\t{}\t{}".format(freq,
                                                                                                     clusters[cluster1][
                                                                                                         0],
                                                                                                     clusters[cluster1][
                                                                                                         1], name,
                                                                                                     gene_pos[name][0],
                                                                                                     gene_pos[name][1],
                                                                                                     clusters[cluster2][
                                                                                                         0],
                                                                                                     clusters[cluster2][
                                                                                                         1], "*","*","*"), file=fout)

                else:

                    for interval in res2:
                        name = interval.data

                        print("t\t{}\talign1\t{}\t{}\t{}\t{}\t{}\talign2\t{}\t{}\t{}\t{}\t{}".format(freq,
                                                                                                     clusters[cluster2][
                                                                                                         0],
                                                                                                     clusters[cluster2][
                                                                                                         1], name,
                                                                                                     gene_pos[name][0],
                                                                                                     gene_pos[name][1],
                                                                                                     clusters[cluster1][
                                                                                                         0],
                                                                                                     clusters[cluster1][
                                                                                                         1], "*", "*",
                                                                                                     "*"), file=fout)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of input file", type=str)
    parser.add_argument("annotation", help="path of file that has annotation", type=str)
    parser.add_argument("output", help="path of file that has annotation", type=str)
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

    cl_tree = intervaltree.IntervalTree()
    for i, cluster in enumerate(clusters):
        cl_tree[cluster[0]:cluster[1]] = i

    count = np.zeros((len(clusters), len(clusters)))
    for aligns in pos_aligns:
        x = cl_tree[aligns[1]:aligns[1] + 150]
        y = cl_tree[aligns[2]:aligns[2] + 150]

        if (len(x) == 0 and len(y) == 0):
            continue

        for x0 in x:
            for y0 in y:
                count[x0.data][y0.data] += 1

    count_result = []
    for ix, iy in np.ndindex(count.shape):
        if count[ix][iy] != 0:
            count_result.append((count[ix][iy], ix, iy))

    count_result.sort(reverse=True)

    gene_tree, gene_pos = build_intervaltree(args.annotation)

    check_annotation(gene_tree, gene_pos, count_result, clusters, args.output)