import argparse
import sys
import os
import pickle


def load_obj(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)


def separate_TP_FP(input_file, dict):
    tp = []
    fp = []
    with open(input_file) as f:
        for line in f:
            line = line.rstrip()
            line_sp = line.split("\t")
            query_id = line_sp[0]
            target_id = line_sp[1]
            if dict.get((query_id, target_id)) or dict.get((target_id, query_id)):
                tp.append(line)

            else:
                fp.append(line)

    return tp, fp


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path of input file", type=str)
    parser.add_argument("dict", help="path of dict file", type=str)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    true_dict = load_obj(args.dict)
    tp, fp = separate_TP_FP(args.input, true_dict)

    save_dir = os.path.dirname(args.input)
    file_name = os.path.basename(args.input)[:-3]
    tp_path = os.path.join(save_dir, file_name + "_tp.out")
    fp_path = os.path.join(save_dir, file_name + "_fp.out")

    with open(tp_path, "w") as fout:
        for ln in tp:
            print(ln, file=fout)

    with open(fp_path, "w") as fout:
        for ln in fp:
            print(ln, file=fout)
