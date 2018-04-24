import argparse
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of input file", type=str)
    parser.add_argument("compare", help="path of compared input file", type=str)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    ref_gene_pos = {}
    with open(args.compare) as f:
        for line in f:
            line = line.strip()
            sp = line.split()
            ref_gene_pos[sp[5]] = (int(sp[9]), int(sp[10]))


    with open(args.input) as f:
        for line in f:
            line = line.strip()
            sp = line.split()
            gene = sp[5]
            x = range(int(sp[9]), int(sp[10]))
            y = range(ref_gene_pos[gene][0], ref_gene_pos[gene][1])
            ovelap_len = len(set(x) & set(y))
            if ovelap_len == 0:
                print(line)