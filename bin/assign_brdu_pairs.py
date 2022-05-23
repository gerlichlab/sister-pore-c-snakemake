import argparse
import pickle
import pandas as pd
import pairtools

def assign_chunk(df, label_lib):
    index_1 = df.read_name.str.cat([df.align1_chrom, df.align1_start.astype("str"), df.align1_end.astype("str")], sep="_")
    index_2 = df.read_name.str.cat([df.align2_chrom, df.align2_start.astype("str"), df.align2_end.astype("str")], sep="_")
    label_1 = index_1.apply(lambda index: label_lib.get(index, None))
    label_2 = index_2.apply(lambda index: label_lib.get(index, None))
    return pd.DataFrame(
        {
            "label_1": label_1,
            "label_2": label_2
        }
    )


# parse args

parser = argparse.ArgumentParser(description='Assign a read to be labelled or unlabelled in the output pairs file depending on the labelling state')
parser.add_argument('--contacts')
parser.add_argument('--pairs')
parser.add_argument('--label_lib')
parser.add_argument('--output')
args = parser.parse_args()

# load label lib

with open(args.label_lib, 'rb') as handle:
    label_library = pickle.load(handle)


# open parquet

contacts = pd.read_parquet(args.contacts)

# assign chunks

result = assign_chunk(contacts, label_library)

# open pairs file to get header

header, pairs_body = pairtools._headerops.get_header(
    pairtools._fileio.auto_open(args.pairs, "r")
)


# extract column names from header

cols = pairtools._headerops.extract_column_names(header)

# construct output pairs

contacts.loc[:, "align1_pos"] = (contacts.align1_start + contacts.align1_end)//2
contacts.loc[:, "align2_pos"] = (contacts.align2_start + contacts.align2_end)//2
contacts.loc[:, "align1_size"] = contacts.align1_end - contacts.align1_start
contacts.loc[:, "align2_size"] = contacts.align2_end - contacts.align2_start

pair_frame = contacts[["read_name", "align1_chrom", "align1_pos", "align2_chrom", "align2_pos", "align1_strand", "align2_strand", "align1_align_score", "align2_align_score", "align1_size", "align2_size"]]
pair_frame.loc[:, "label_1"] = result["label_1"]
pair_frame.loc[:, "label_2"] = result["label_2"]

pair_frame = pair_frame.rename(columns={"read_name": "readID", "align1_chrom": "chr1", "align1_pos": "pos1",  "align2_chrom": "chr2", "align2_pos": "pos2", "align1_strand": "strand1", "align2_strand": "strand2"})

pair_frame_clean = pair_frame.dropna()

new_cols = f'#columns: {" ".join(list(pair_frame_clean.columns))}'

header[-1] = new_cols

# write output

out_stream = pairtools._fileio.auto_open(args.output, nproc=3, mode="w")

out_stream.writelines([line + "\n" for line in header])

pair_frame_clean.to_csv(out_stream, sep="\t", columns=list(pair_frame_clean.columns), header=False, index=False)


