import argparse
import sys
from collections import defaultdict

import intervaltree
import numpy as np
import pandas as pd

import SAMparser


def build_intervaltree(input, input2=None, rc=False):
    gene_pos = {}
    tree = intervaltree.IntervalTree()
    with open(input) as f:
        for line in f:
            if line[0] != "#":
                sp = line.strip().split()
                if sp[3] == "277694":
                    if (sp[6] == '-' and rc) or (sp[6] == '+' and not rc):
                        start = min(int(sp[5]), int(sp[4]))
                        end = max(int(sp[5]), int(sp[4]))
                        name = sp[8]
                        tree[start:end + 1] = name
                        gene_pos[name] = (start, end + 1)

    if input2:
        with open(input2) as f:
            for line in f:
                sp = line.strip().split()
                if (sp[5] == '-' and rc) or (sp[5] == '+' and not rc):
                    start = min(int(sp[1]), int(sp[2]))
                    end = max(int(sp[1]), int(sp[2]))
                    name = sp[3]
                    tree[start:end + 1] = name
                    gene_pos[name] = (start, end + 1)

    return tree, gene_pos


def check_annotation(tree, gene_pos, count, clusters):
    label = []
    freqs = []
    clus1_start = []
    clus1_end = []
    gene1 = []
    gene1_start = []
    gene1_end = []
    clus2_start = []
    clus2_end = []
    gene2 = []
    gene2_start = []
    gene2_end = []

    for ct in count:

        cluster1 = ct[1]
        cluster2 = ct[2]
        freq = ct[0]/2
        res1 = tree[clusters[cluster1][0]:clusters[cluster1][1]]
        res2 = tree[clusters[cluster2][0]:clusters[cluster2][1]]


        if (len(res1) == 0 and len(res2) == 0):
            freqs.append(freq)
            label.append('f')
            clus1_start.append(clusters[cluster1][0])
            clus1_end.append(clusters[cluster1][1])
            gene1.append(None)
            gene1_start.append(None)
            gene1_end.append(None)
            clus2_start.append(clusters[cluster2][0])
            clus2_end.append(clusters[cluster2][1])
            gene2.append(None)
            gene2_start.append(None)
            gene2_end.append(None)


        elif (len(res1) != 0 and len(res2) != 0):
            for interval1 in res1:
                for interval2 in res2:
                    name1 = interval1.data
                    name2 = interval2.data
                    if 'nsRNA' in name1:
                        name1, name2 = name2, name1
                        clu1, clu2 = cluster2, cluster1
                        freqs.append(freq)
                        label.append('t')
                        clus1_start.append(clusters[clu1][0])
                        clus1_end.append(clusters[clu1][1])
                        gene1.append(name1)
                        gene1_start.append(gene_pos[name1][0])
                        gene1_end.append(gene_pos[name1][1])
                        clus2_start.append(clusters[clu2][0])
                        clus2_end.append(clusters[clu2][1])
                        gene2.append(name2)
                        gene2_start.append(gene_pos[name2][0])
                        gene2_end.append(gene_pos[name2][1])

                    elif 'nsRNA' in name2:
                        freqs.append(freq)
                        label.append('t')
                        clus1_start.append(clusters[cluster1][0])
                        clus1_end.append(clusters[cluster1][1])
                        gene1.append(name1)
                        gene1_start.append(gene_pos[name1][0])
                        gene1_end.append(gene_pos[name1][1])
                        clus2_start.append(clusters[cluster2][0])
                        clus2_end.append(clusters[cluster2][1])
                        gene2.append(name2)
                        gene2_start.append(gene_pos[name2][0])
                        gene2_end.append(gene_pos[name2][1])


                    else:
                        freqs.append(freq)
                        label.append('f')
                        clus1_start.append(clusters[cluster1][0])
                        clus1_end.append(clusters[cluster1][1])
                        gene1.append(name1)
                        gene1_start.append(gene_pos[name1][0])
                        gene1_end.append(gene_pos[name1][1])
                        clus2_start.append(clusters[cluster2][0])
                        clus2_end.append(clusters[cluster2][1])
                        gene2.append(name2)
                        gene2_start.append(gene_pos[name2][0])
                        gene2_end.append(gene_pos[name2][1])
        else:
            if len(res1) != 0:
                for interval in res1:
                    name = interval.data

                    freqs.append(freq)
                    label.append('u')
                    clus1_start.append(clusters[cluster1][0])
                    clus1_end.append(clusters[cluster1][1])
                    gene1.append(name)
                    gene1_start.append(gene_pos[name][0])
                    gene1_end.append(gene_pos[name][1])
                    clus2_start.append(clusters[cluster2][0])
                    clus2_end.append(clusters[cluster2][1])
                    gene2.append(None)
                    gene2_start.append(None)
                    gene2_end.append(None)

            else:

                for interval in res2:
                    name = interval.data

                    freqs.append(freq)
                    label.append('u')
                    clus1_start.append(clusters[cluster2][0])
                    clus1_end.append(clusters[cluster2][1])
                    gene1.append(name)
                    gene1_start.append(gene_pos[name][0])
                    gene1_end.append(gene_pos[name][1])
                    clus2_start.append(clusters[cluster1][0])
                    clus2_end.append(clusters[cluster1][1])
                    gene2.append(None)
                    gene2_start.append(None)
                    gene2_end.append(None)

    df = pd.DataFrame({'label': label, 'freq': freqs, 'cluster1_start': clus1_start, 'cluster1_end': clus1_end,
                       'gene1': gene1, 'gene1_start': gene1_start, 'gene1_end': gene1_end,
                       'cluster2_start': clus2_start, 'cluster2_end': clus2_end,
                       'gene2': gene2, 'gene2_start': gene2_start, 'gene2_end': gene2_end})

    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of input file", type=str)
    parser.add_argument("annotation", help="path of file that has annotation", type=str)
    parser.add_argument("output", help="path of output", type=str)
    parser.add_argument("--reverse", help="strand", default=False, action='store_true')
    parser.add_argument("--annotation2", help="path of small RNA file", default=None, type=str)
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

                if args.reverse == record.rc:

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
                    if record.sa_rname == record.rname and args.reverse == record.sa_strand:
                        pos.append(record.sa_pos)

        min_pos = np.min(pos)
        max_pos = np.max(pos)

        if max_pos - min_pos >= 1000:
            pos_aligns.append((min_pos, min_pos, max_pos))
            pos_aligns.append((max_pos, min_pos, max_pos))

    clusters = [[pos_aligns[0][0], pos_aligns[0][0] + 150]]
    last = 0
    pos_aligns.sort()
    for aligns in pos_aligns:
        align_pos, _, _ = aligns
        if align_pos + 150 >= last:
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

    gene_tree, gene_pos = build_intervaltree(args.annotation, args.annotation2, args.reverse)

    result = check_annotation(gene_tree, gene_pos, count_result, clusters)

    result = result.sort_values(['gene2', 'freq'], ascending=[True, False])
    result.to_csv(args.output, columns=['label', 'freq', 'cluster1_start', 'cluster1_end',
                       'gene1', 'gene1_start', 'gene1_end',
                       'cluster2_start', 'cluster2_end',
                       'gene2', 'gene2_start', 'gene2_end'], sep='\t', index=False)
