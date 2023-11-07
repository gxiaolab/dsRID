import sys
sys.path.insert(1, '/u/home/r/ryo10244/Xiao_lab/dsRNA_pred/script/utils')

import re
import pysam
import argparse as ap
import pandas as pd
import numpy as np
from utils import *
#import multiprocessing as mp
#from functools import partial

def main(args):
    argp = ap.ArgumentParser(description="extract features for dsrna prediction for real dsrna region",
                             formatter_class=ap.ArgumentDefaultsHelpFormatter)
    argp.add_argument(
        "-b", "--bam_file",
        help="input bam file, with cs tags, sorted and indexed",
        type=str,
        default=""
    )

    argp.add_argument(
        "-o", "--out_file",
        help="output file location",
        type=str,
        default="./data/train.tsv"
    )

    argp.add_argument(
        "-d", "--ds_file",
        help="input double stranded RNA file, tab delimnated bed",
        type=str,
        default="./data/hg38_dsrna.bed"
    )

    argp.add_argument(
        "-r", "--read_threshold",
        help="number of coverage required for prediction",
        type=int,
        default=6
    )

    args = argp.parse_args(args)

    mapq_thr = 20
    #chr = args.chr

    dsrna = pd.read_csv(args.ds_file, sep='\t', names=['chr', 'start', 'end', 'name'], header=None)
    print(dsrna)
    dsrna = dsrna.sort_values(by=['chr', 'start'])
    sam = pysam.AlignmentFile(args.bam_file, 'rb')

    feat_lst = list()

    i = 0
    n_row = dsrna.shape[0]

    for chr, start, end, name in zip(dsrna['chr'], dsrna['start'],
     dsrna['end'], dsrna['name']):
        read_list = list()
        #print(chr, start, end)
        coverage = sam.count(chr, start, end)
        i += 1
        if  coverage < args.read_threshold:
            continue
        for read in sam.fetch(chr, start, end):
            if read.mapq < mapq_thr:
                continue
            elif read.is_secondary: # skip secondary reads
                continue
            # 0 based
            pos = read.reference_start
            read_start = read.reference_start
            read_end = read.reference_end
            cs = cs_to_df(read.get_tag('cs'), pos)
            #print(cs)
            cs_splice = cs.loc[cs['ope'] == '~']
            #print(cs_splice)
            #continue
            for ri, row in cs_splice.iterrows():
                low = int(row['low'])
                high = int(row['high'])
                # if low < start or high > end:
                #     continue
                a = row['val'][0:2]
                b = row['val'][-2:]
                length = int(row['val'][2:-2])
                skipped_bases = read.query_sequence[low - start : high - start]
                read_list.append(
                    [read.query_name,
                    chr, low, high, length,
                    a, b, read_start, read_end, skipped_bases]
                )
        mapped_splices = pd.DataFrame(read_list,
        columns=["read", "chr", "pos_start",
         "pos_end", "pos_len", "bp_start", "bp_end",
         "read_start", "read_end", "skipped_bases"])

        sp_sites = mapped_splices.loc[(mapped_splices['pos_start'] > start - 10) &
      (mapped_splices['pos_end'] < end + 10) ]
        if len(sp_sites) == 0:
            continue
        #print(mapped_splices)
        feat_mat = reads_to_feature(mapped_splices, start, end, coverage)
        feat_mat['name'] = name
        feat_mat['chr'] = chr
        feat_mat['start'] = start
        feat_mat['end'] = end
        feat_lst.append(feat_mat)
        if len(feat_lst) % 1000 == 0:
            feat_total = pd.concat(feat_lst)
            feat_total.to_csv(args.out_file, sep='\t')
        print(len(feat_lst), i, n_row)
    feat_total = pd.concat(feat_lst)
    feat_total.to_csv(args.out_file, sep='\t')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))