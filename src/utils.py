import re
from statistics import mean
import pysam
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

def split_cs_string(cs_string):
    return list(
        zip(
            re.sub('[0-9a-z]', ' ', cs_string).split(),
            re.sub('[:*\-+~]', ' ', cs_string).split()
        )
    )

cslenfuncs = {
    ':': int,
    '*': lambda x: 1,
    '+': lambda x: 0,
    '-': len,
    '~': lambda x: int(re.sub('[a-z]', '', x))
}

cigarlenfuncs = {
    0 : lambda x: x,
    1 : lambda x: 0,
    2 : lambda x: x,
    3 : lambda x: x,
    4 : lambda x: x,
    5 : lambda x: x,
    6 : lambda x: x,
    7 : lambda x: x,
    8 : lambda x: x
}

def cs_to_df(cs_string, pos):
    cs = split_cs_string(cs_string)
    cslist = list()
    for a, b in cs:
        low = pos
        pos += cslenfuncs[a](b)
        high = pos
        cslist.append([low, high, a, b])
    csdf = pd.DataFrame(
        np.row_stack(cslist),
        columns=['low', 'high', 'ope', 'val']
    )
    csdf.loc[:, 'low'] = csdf['low'].astype(int)
    csdf.loc[:, 'high'] = csdf['high'].astype(int)
    return csdf

def CIGAR_to_df(cigar_tuples, pos):
    cigar_list = list()
    for symbol, length in cigar_tuples:
        low = pos
        pos += cigarlenfuncs[symbol](length)
        high = pos
        cigar_list.append([low, high, symbol, length])
    cigar_df = pd.DataFrame(
        np.row_stack(cigar_list),
        columns=['low', 'high', 'ope', 'val']
    )
    return cigar_df

def reads_to_feature(mapped_splices, start, end, coverage):
    sp_sites = mapped_splices.loc[(mapped_splices['pos_start'] > start - 10) &
      (mapped_splices['pos_end'] < end + 10) ]

    group_ids, group_num, group_std = get_sp_cluster(sp_sites)

    mean_sp_start = np.mean(sp_sites['pos_start'])
    mean_sp_end = np.mean(sp_sites['pos_end'])
    var_sp_start = np.std(sp_sites['pos_start'] - mean_sp_start)
    var_sp_end = np.std(sp_sites['pos_end'] - mean_sp_end)

    res = pd.DataFrame(data ={"std_start" : var_sp_start,
     "std_end": var_sp_end, "mean_start":mean_sp_start,
     "mean_end": mean_sp_end, "len_skip":mean_sp_end - mean_sp_start,
      "coverage" : coverage, "num_skip": sp_sites.shape[0],
      "skip_ratio": sp_sites.shape[0] / coverage,
       "group_num": group_num, "group_std" : group_std 
       }, index=[0], 
       )
    sp_sites['group'] = group_ids

    #print(sp_sites['skipped_bases'])

    gc_skip = np.mean(sp_sites["skipped_bases"].apply(lambda st : len(re.findall("[GC]", st)) / len(st) if len(st) > 0 else 0 ))

    res['gc_skip'] = gc_skip

    count_bp_start = count_bp(sp_sites['bp_start'], True)
    count_bp_end = count_bp(sp_sites['bp_end'], False)

    res = pd.concat([res, count_bp_start, count_bp_end], axis=1)

    return res

def get_sp_cluster(sp_sites, window_size=50):
    "Given splicing sites return groups of sp sites that are at least 100 bps apart"
    group_id = 0
    pos_temp = min(sp_sites['pos_start'])
    group_ids = [0]
    group_dic = {group_id : [pos_temp]}
    for pos_start in sp_sites['pos_start'][1:].sort_values():
        if pos_start - pos_temp > window_size:
            group_id += 1
            group_dic[group_id] =  [pos_start]
            pos_temp = mean(group_dic[group_id])
        else:
            group_dic[group_id].append(pos_start)
            pos_temp = mean(group_dic[group_id])
        group_ids.append(group_id)
    
    stds = []
    for key in group_dic.keys():
        std_temp = np.std(group_dic[key])
        stds.append(std_temp)
    return group_ids, max(group_ids) + 1, mean(stds)


def count_bp(bps, start):
    res = dict()
    prefix = 'bp_start_' if start else 'bp_end_'
    for bp in bps.unique():
        res[prefix + bp] = 0
    for bp in bps:
        res[prefix + bp] += 1 / len(bps)
    return pd.DataFrame(data=res, index=[0])

def get_gtf_splice_pos(gtf_file, chromosomes):
    gtf = pysam.TabixFile(gtf_file)

    poslist = defaultdict(list)

    for chrom in chromosomes:
        for gtf_entry in gtf.fetch(
                chrom,
                parser=pysam.asGTF()):
            if gtf_entry.feature == 'exon':
                poslist[chrom].append((gtf_entry.start, gtf_entry.end))
    for key in poslist.keys():
        poslist[key].sort()
    return poslist
