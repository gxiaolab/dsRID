import sys
sys.path.insert(1, '/u/home/r/ryo10244/Xiao_lab/dsRNA_pred/script/utils')

import re
import pysam
import argparse as ap
import pandas as pd
import numpy as np
from utils import *
import random
from os.path import exists
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
        default="/data/dsRID_whole.tsv"
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
    window = 2500
    slide = 1250
    #chrs = args.chr
    chrs = ["chrX", "chrY"]
    chrs.extend(["chr" + str(i) for i in range(1, 23)])
    print(chrs)
    sample_num = 1000000
    chr_max = {}
    for chr in chrs:
        chr_max[chr] = 0

    sam = pysam.AlignmentFile(args.bam_file, 'rb')
    feat_lst = list()

    prev_exists = exists(args.out_file)

    frame_prev = pd.DataFrame()

    if prev_exists:
        frame_prev = pd.read_csv(args.out_file, sep='\t')
        for chr in chrs:
            frame_prev_chr = frame_prev.loc[frame_prev['chr'] == chr]
            if frame_prev_chr.shape[0] > 0:
                chr_max[chr] = max(frame_prev_chr['end'])
    #sp_dic = get_gtf_splice_pos(args.splice_anno, chrs)
    i = 0
    start = 0
    end = start + window
    for chr in chrs:
        start = 0
        end = start + window
        name = chr + ":" + str(start) + "-" + str(end)
        
        while(end < 2.5e8):
            if end < chr_max[chr]:
                end = chr_max[chr]
                start = end - window
                name = chr + ":" + str(start) + "-" + str(end)
                continue
            #print(chr, start, end)
            
            coverage = sam.count(chr, start, end)
            read_list = list()
            
            if  coverage < args.read_threshold:
                start = start + slide
                end = start + window
                name = chr + ":" + str(start) + "-" + str(end)
                continue
            for read in sam.fetch(chr, start, end, until_eof=True):
                if read.mapq < mapq_thr:
                    continue
                elif read.is_secondary: # skip secondary reads
                    continue
                # 0 based
                pos = read.reference_start
                read_start = read.reference_start
                read_end = read.reference_end
                # print(read.cigarstring)
                # print(read.cigartuples)
                # print(read.get_tags())
                cs = CIGAR_to_df(read.cigartuples, pos)
                cs_splice = cs.loc[cs['ope'] == 3]
                #print(cs_splice)
                #continue
                for ri, row in cs_splice.iterrows():
                    low = int(row['low'])
                    high = int(row['high'])
                    # if low < start or high > end:
                    #     continue
                    bp_start = read.query_sequence[low-read_start - 2 : low - read_start]
                    bp_end = read.query_sequence[high-read_start : high - read_start + 2]
                    #print(bp_start, bp_end)
                    length = int(row['val'])
                    skipped_bases = read.query_sequence[low - read_start : high - read_start]
                    read_list.append(
                        [read.query_name,
                        chr, low, high, length,
                        bp_start, bp_end, read_start, read_end, skipped_bases]
                    )
            mapped_splices = pd.DataFrame(read_list,
            columns=["read", "chr", "pos_start",
            "pos_end", "pos_len", "bp_start", "bp_end",
            "read_start", "read_end", "skipped_bases"])

            sp_sites = mapped_splices.loc[(mapped_splices['pos_start'] > start - 10) &
        (mapped_splices['pos_end'] < end + 10) ]
            if len(sp_sites) == 0:
                start = start + slide
                end = start + window
                continue
            #print(mapped_splices)
            feat_mat = reads_to_feature(mapped_splices, start, end, coverage)
            feat_mat['name'] = name
            feat_mat['chr'] = chr
            feat_mat['start'] = start
            feat_mat['end'] = end
            feat_lst.append(feat_mat)
            if len(feat_lst) % 1000 == 0:
                if prev_exists:
                    feat_total = pd.concat(feat_lst)
                    feat_total = pd.concat([feat_total, frame_prev])
                    feat_total.to_csv(args.out_file, sep='\t')
                else:
                    feat_total = pd.concat(feat_lst)
                    feat_total.to_csv(args.out_file, sep='\t')
            #print(len(feat_lst))
            print(i, "/", sample_num, name)
            i += 1
            start = end
            end = start + window
            name = chr + ":" + str(start) + "-" + str(end)
    feat_total = pd.concat(feat_lst)
    feat_total = pd.concat([feat_total, frame_prev])
    feat_total.to_csv(args.out_file, sep='\t')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))