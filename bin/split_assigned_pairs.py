import argparse
from pathlib import Path
import pandas as pd
import pairtools


def write_pairs(header, body, outpath):
    out_stream = pairtools._fileio.auto_open(outpath, nproc=3, mode="w")
    out_stream.writelines([line + "\n" for line in header])
    body.to_csv(
        out_stream, sep="\t", columns=list(body.columns), header=False, index=False
    )


def split_pairs(pairs):
    trusty_pairs = pairs.loc[
        (pairs.align1_size > 100) & (pairs.align2_size > 100), :
    ].rename(columns={"chr1": "chrom1", "chr2": "chrom2"})
    # define constants
    IS_TRANS = (
        (
            (trusty_pairs.label_1 & trusty_pairs.strand1)
            & (trusty_pairs.label_2 & ~trusty_pairs.strand2)
        )
        | (
            (trusty_pairs.label_1 & trusty_pairs.strand1)
            & (~trusty_pairs.label_2 & trusty_pairs.strand2)
        )
        | (
            (~trusty_pairs.label_1 & ~trusty_pairs.strand1)
            & (~trusty_pairs.label_2 & trusty_pairs.strand2)
        )
        | ((~trusty_pairs.label_1 & ~trusty_pairs.strand1))
        & (trusty_pairs.label_2 & ~trusty_pairs.strand2)
        | (
            (trusty_pairs.label_1 & ~trusty_pairs.strand1)
            & (trusty_pairs.label_2 & trusty_pairs.strand2)
        )
        | (
            (~trusty_pairs.label_1 & trusty_pairs.strand1)
            & (trusty_pairs.label_2 & trusty_pairs.strand2)
        )
        | (
            (~trusty_pairs.label_1 & trusty_pairs.strand1)
            & (~trusty_pairs.label_2 & ~trusty_pairs.strand2)
        )
        | (
            (trusty_pairs.label_1 & ~trusty_pairs.strand1)
            & (~trusty_pairs.label_2 & ~trusty_pairs.strand2)
        )
    )
    # get cis and trans contacts
    trans_pairs = trusty_pairs.loc[IS_TRANS, :]
    cis_pairs = trusty_pairs.loc[~IS_TRANS, :]
    return trusty_pairs, cis_pairs, trans_pairs


# parse args

parser = argparse.ArgumentParser(
    description="Split assigned pairs into different types."
)
parser.add_argument("--pairs")
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

all_pairs, cis_pairs, trans_pairs = split_pairs(pairs)

write_pairs(header, all_pairs, base_name + "all_reads.cis_and_trans.pairs.gz")
write_pairs(header, cis_pairs, base_name + "all_reads.cis.pairs.gz")
write_pairs(header, trans_pairs, base_name + "all_reads.trans.pairs.gz")

# generate output without two unlabelled contacts

all_wo_unlabelled, cis_wo_unlabelled, trans_wo_unlabelled = split_pairs(
    pairs.loc[~(~pairs.label_1 & ~pairs.label_2), :]
)
write_pairs(header, all_wo_unlabelled, base_name + "no_unlabelled_c.cis_and_trans.pairs.gz")
write_pairs(
    header, cis_wo_unlabelled, base_name + "no_unlabelled_c.cis.pairs.gz"
)
write_pairs(header, trans_wo_unlabelled, base_name + "no_unlabelled_c.trans.pairs.gz")

# generate output with mapping quality heuristic filter

all_pairs_mq, cis_mq, trans_mq = split_pairs(
    pairs.loc[
        ((pairs.align1_mapping_quality < 40) | (pairs.align1_mapping_quality > 55))
        & ((pairs.align2_mapping_quality < 40) | (pairs.align2_mapping_quality > 55)),
        :
    ]
)
write_pairs(header, all_pairs_mq, base_name + "mq_heuristic.cis_and_trans.pairs.gz")
write_pairs(header, cis_mq, base_name + "mq_heuristic.cis.pairs.gz")
write_pairs(header, trans_mq, base_name + "mq_heuristic.trans.pairs.gz")
