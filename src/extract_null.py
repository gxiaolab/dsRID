import sys

from scipy.fft import skip_backend
sys.path.insert(1, '/u/home/r/ryo10244/Xiao_lab/dsRNA_pred/script/utils')

import re
import pysam
import argparse as ap
import pandas as pd
import numpy as np
from utils import *
import random
#import multiprocessing as mp
#from functools import partial

def main(args):
    argp = ap.ArgumentParser(description="extract features for dsrna prediction from randomly sampled region",
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
        default="./data/null.tsv"
    )

    argp.add_argument(
        "-r", "--read_threshold",
        help="number of reads required for prediction",
        type=int,
        default=6
    )

    argp.add_argument(
        "-s", "--splice_anno",
        help="gtf file for splicing region",
        type=str,
        default="./data/gencode.v34.annotation.sorted.gtf.gz"
    )

    argp.add_argument(
        "-c", "--chr",
        help = "chromosomes to be analyzed",
        nargs='*',
        type=str,
        default = 'chr1'
    )

    args = argp.parse_args(args)

    mapq_thr = 20
    window = 2500
    #chrs = args.chr
    chrs = ["chrX", "chrY"]
    chrs.extend(["chr" + str(i) for i in range(1, 23)])
    print(chrs)
    sample_num = 10000

    sam = pysam.AlignmentFile(args.bam_file, 'rb')
    feat_lst = list()

    sp_dic = get_gtf_splice_pos(args.splice_anno, chrs)

    while (len(feat_lst) < sample_num):
        chr = random.choice(chrs)
        exonsite_ind = random.choice(range(len(sp_dic[chr])))
        try:
            int_start = sp_dic[chr][exonsite_ind][1]
            int_end = sp_dic[chr][exonsite_ind + 1][0]
        except:
            continue
        if int_start + 10 - window > int_end - 10 + window:
            continue
        start = np.random.randint(int_start + 10 - window, int_end - 10 + window)
        end = start + window
        print(chr, start, end)
        coverage = sam.count(chr, start, end)
        read_lst = list()
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
            cs_splice = cs.loc[cs['ope'] == '~']
            #print(cs_splice)
            if len(cs_splice) == 0:
                continue
            for ri, row in cs_splice.iterrows():
                low = int(row['low'])
                high = int(row['high'])
                # if low < start or high > end:
                #     continue
                a = row['val'][0:2]
                b = row['val'][-2:]
                length = int(row['val'][2:-2])
                skipped_bases = read.query_sequence[low - start : high - start]
                #print(read.query_sequence)
                #print(skipped_bases)
                read_lst.append(
                    [read.query_name,
                    chr, low, high, length,
                    a, b, read_start, read_end, skipped_bases]
                )
        
        mapped_splices = pd.DataFrame(read_lst,
        columns=["read", "chr", "pos_start",
        "pos_end", "pos_len", "bp_start", "bp_end",
        "read_start", "read_end", "skipped_bases"])

        sp_sites = mapped_splices.loc[(mapped_splices['pos_start'] > start - 10) &
      (mapped_splices['pos_end'] < end + 10) ]
        if len(sp_sites) == 0:
            continue

        print(sp_sites['pos_start'][sp_sites["pos_start"].isna()])
        feat_mat = reads_to_feature(mapped_splices, start, end, coverage)
        feat_mat['name'] = chr + ":" + str(start) + "-" + str(end)
        feat_mat['chr'] = chr
        feat_mat['start'] = start
        feat_mat['end'] = end
        print(feat_mat)
        feat_lst.append(feat_mat)
        if len(feat_lst) % 1000 == 0:
            feat_total = pd.concat(feat_lst)
            feat_total.to_csv(args.out_file, sep='\t')
        print(len(feat_lst))
    feat_total = pd.concat(feat_lst)
    feat_total.to_csv(args.out_file, sep='\t')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))