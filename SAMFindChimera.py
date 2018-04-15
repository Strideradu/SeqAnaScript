import argparse
import sys
import os
from collections import defaultdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of input file", type=str)
    parser.add_argument("output", help="path of output file", type=str)
    parser.add_argument("id_list", help="path of id list file", type=str)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    aligns = defaultdict(list)
    with open(args.input) as f:
        for line in f:
            if line[0] != "@":
                line = line.rstrip()
                sp = line.split('\t')
                if sp[2]!='*':
                    aligns[sp[0]].append(line)
    with open(args.output,'w') as f1:
        with open(args.id_list, 'w') as f2:
            for key, value in aligns.items():
                if len(value) > 2:
                    for line in value:
                        print(line, file=f1)

                    print(key, file=f2)


