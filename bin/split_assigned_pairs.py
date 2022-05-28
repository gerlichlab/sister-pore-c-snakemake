import argparse
from pathlib import Path
import pandas as pd
import pairtools

def write_pairs(header, body, outpath):
    out_stream = pairtools._fileio.auto_open(outpath, nproc=3, mode="w")
    out_stream.writelines([line + "\n" for line in header])
    body.to_csv(out_stream, sep="\t", columns=list(body.columns), header=False, index=False)

# parse args

parser = argparse.ArgumentParser(description='Split assigned pairs into different types.')
parser.add_argument('--pairs')
args = parser.parse_args()

base_name = args.pairs.split("brdu_assigned.pairs.gz")[0]

# read in pairs

header, pairs_body = pairtools._headerops.get_header(
    pairtools._fileio.auto_open(args.pairs, "r")
)

cols = pairtools._headerops.extract_column_names(header)

pairs = pd.read_csv(pairs_body, sep="\t", names=cols)

# generate labelled_only output

labelled_pairs = pairs.loc[(pairs.label_1 & pairs.label_2), :]
trans_pairs = labelled_pairs.loc[(labelled_pairs.strand1 != labelled_pairs.strand2), :]
cis_pairs = labelled_pairs.loc[(labelled_pairs.strand1 == labelled_pairs.strand2), :]

write_pairs(header, labelled_pairs, base_name + "labelled_only.cis_and_trans.pairs.gz")
write_pairs(header, cis_pairs, base_name + "labelled_only.cis.pairs.gz")
write_pairs(header, trans_pairs, base_name + "labelled_only.trans.pairs.gz")

# generate all output

trusty_pairs = pairs.loc[(pairs.align1_size > 100) & (pairs.align2_size > 100), :]

is_trans = (
    ((trusty_pairs.label_1 & trusty_pairs.strand1) & (trusty_pairs.label_2 & ~trusty_pairs.strand2)) |
    ((trusty_pairs.label_1 & trusty_pairs.strand1) & (~trusty_pairs.label_2 & trusty_pairs.strand2)) |
    ((~trusty_pairs.label_1 & ~trusty_pairs.strand1) & (~trusty_pairs.label_2 & trusty_pairs.strand2)) |
    ((~trusty_pairs.label_1 & ~trusty_pairs.strand1)) & (trusty_pairs.label_2 & ~trusty_pairs.strand2) |
    ((trusty_pairs.label_1 & ~trusty_pairs.strand1) & (trusty_pairs.label_2 & trusty_pairs.strand2)) |
    ((~trusty_pairs.label_1 & trusty_pairs.strand1) & (trusty_pairs.label_2 & trusty_pairs.strand2)) |
    ((~trusty_pairs.label_1 & trusty_pairs.strand1) & (~trusty_pairs.label_2 & ~trusty_pairs.strand2)) |
    ((trusty_pairs.label_1 & ~trusty_pairs.strand1) & (~trusty_pairs.label_2 & ~trusty_pairs.strand2)) 
)

trans_pairs = trusty_pairs.loc[is_trans, :]
cis_pairs = trusty_pairs.loc[~is_trans, :]

write_pairs(header, trusty_pairs, base_name + "all_reads.cis_and_trans.pairs.gz")
write_pairs(header, cis_pairs, base_name + "all_reads.cis.pairs.gz")
write_pairs(header, trans_pairs, base_name + "all_reads.trans.pairs.gz")