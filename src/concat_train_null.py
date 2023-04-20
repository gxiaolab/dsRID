import sys
sys.path.insert(1, '/u/home/r/ryo10244/Xiao_lab/dsRNA_pred/script/utils')

import re
import pysam
import argparse as ap
import pandas as pd
import numpy as np
from utils import *

def main(args):
    argp = ap.ArgumentParser(description="concatanate train and null file",
                             formatter_class=ap.ArgumentDefaultsHelpFormatter)
    argp.add_argument(
        "-t", "--train",
        help="input train file",
        type=str,
        default="/data/dsRID_train.tsv"
    )

    argp.add_argument(
        "-n", "--null",
        help="input null file",
        type=str,
        default="/data/dsRID_null.tsv"
    )

    argp.add_argument(
        "-o", "--out_file",
        help="output file location after concatanation",
        type=str,
        default="/data/dsRID_data.tsvv"
    )
    args = argp.parse_args(args)


    dsrna = pd.read_csv(args.train, sep='\t', index_col=0)

    ds_filter = ~dsrna.columns.isin(["mean_start", "mean_end",
    "coverage", "num_skip", "name",  "chr", "start", "end", "label"])

    dsrna = dsrna.fillna(value=0)

    #dsrna_sub = dsrna.loc[:, ds_filter]

    #dsrna = dsrna.loc[dsrna_sub.sum(axis=1) != 0, :]

    #print(dsrna)

    dsrna['label'] = 1

    null = pd.read_csv(args.null, sep='\t', index_col=0)

    null['label'] = 0

    all = pd.concat([dsrna, null], axis=0)

    all.to_csv(args.out_file,
    sep='\t', index = False)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))