import argparse
import sys
from Bio import AlignIO

def calErrorRates(file_path):
    bases = 0
    errors = 0
    insertions = 0
    deletions = 0
    for multiple_alignment in AlignIO.parse(file_path, "maf"):
        multiple_alignment = list(multiple_alignment)
        bases += multiple_alignment[0].annotations["size"]
        for i in range(len(multiple_alignment[0].seq)):
            if multiple_alignment[0].seq[i]!=multiple_alignment[1].seq[i]:
                errors+=1
                if multiple_alignment[0].seq[i] == '-':
                    insertions += 1

                if multiple_alignment[1].seq[i] == '-':
                    deletions += 1

    print('Error rates is {}'.format(errors/bases))
    print('Insertion rates is {}'.format(insertions/bases))
    print('Deletion rates is {}'.format(deletions / bases))
    print('Substitution rates is {}'.format((errors - insertions -deletions)/bases))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of input file", type=str)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    calErrorRates(args.input)
