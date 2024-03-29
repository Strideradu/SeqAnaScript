import dill as pickle
import argparse
import sys
import intervaltree
from collections import defaultdict


def load_obj(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)


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


def check_annotation(tree, gene_pos, input, dict, output):
    with open(args.output, "w") as fout:
        with open(input) as f:
            for line in f:
                sp = line.strip().split()
                if sp[3] == "True":
                    x = int(sp[0])
                    y = int(sp[1])
                    print(file=fout)
                    print(sp[2],file=fout)
                    result = {}
                    for align in dict[x][y]:
                        res1 = tree[align[0]:align[0] + 150]
                        res2 = tree[align[1]:align[1] + 150]

                        if (len(res1) == 0 and len(res2) == 0) or (len(res1) != 0 and len(res2) != 0):
                            continue

                        else:
                            if len(res1) != 0:
                                for interval in res1:
                                    name = interval.data
                                    if name not in result:
                                        result[name] = (align[1], align[1] + 150)

                                    result[name] = (min(align[1], result[name][0]), max(align[1] + 150, result[name][1]))

                            else:

                                for interval in res2:
                                    name = interval.data
                                    if name not in result:
                                        result[name] = (align[0], align[0] + 150)

                                    result[name] = (min(align[0], result[name][0]), max(align[0] + 150, result[name][1]))
                    for gene, pos in result.items():
                        print("{}\t{}\t{}\t{}\t{}".format(gene, gene_pos[gene][0], gene_pos[gene][1], pos[0], pos[1]), file=fout)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input1", help="path of file that has the number of pairs", type=str)
    parser.add_argument("input2", help="path of file that has align pos dict", type=str)
    parser.add_argument("annotation", help="path of file that has annotation", type=str)
    parser.add_argument("output", help="path of file that has annotation", type=str)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    aligns = load_obj(args.input2)
    tree, gene_pos = build_intervaltree(args.annotation)

    result = check_annotation(tree, gene_pos, args.input1, aligns, args.output)