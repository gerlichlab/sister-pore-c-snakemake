import pickle
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sbn
from pathlib import Path
import pairtools
import datetime
import argparse

# constants

BATCH_ID = r'(batch\d{1,2})'

# parse arguments

parser = argparse.ArgumentParser(description='Create pore-c qc plots')
parser.add_argument('--detect_paths', nargs="+")
parser.add_argument('--label_library_paths', nargs="+")
parser.add_argument('--all_pairs_paths', nargs="+")
parser.add_argument('--all_reads_cis_trans', nargs="+")
parser.add_argument('--all_reads_cis', nargs="+")
parser.add_argument('--all_reads_trans', nargs="+")
parser.add_argument('--labelled_cis_trans', nargs="+")
parser.add_argument('--labelled_cis', nargs="+")
parser.add_argument('--labelled_trans', nargs="+")
parser.add_argument('--no_c_cis_trans', nargs="+")
parser.add_argument('--no_c_cis', nargs="+")
parser.add_argument('--no_c_trans', nargs="+")
parser.add_argument('--mq_cis_trans', nargs="+")
parser.add_argument('--mq_cis', nargs="+")
parser.add_argument('--mq_trans', nargs="+")
parser.add_argument('--output')
args = parser.parse_args()

# Functions

def make_title(title, ax, subtitle=None):
    plt.sca(ax)
    plt.axis("off")
    date = datetime.date.today()
    ax.text(
        0.5,
        0.8,
        f"{title} from {date}",
        fontsize=24,
        horizontalalignment="center",
        verticalalignment="center",
    )
    if subtitle is not None:
        ax.text(
            0,
            0.3,
            subtitle,
            fontsize=18,
            horizontalalignment="left",
            verticalalignment="center",
        )

class Sample():
    def __init__(self, name):
        self.name = name
        self._reads = []
        self._read_names = set()
    def add_read(self, read):
        if read._read_id not in self._read_names:
            self._reads.append(read)
            self._read_names.add(read._read_id)
    def get_thymidines(self):
        output = []
        for read in self._reads:
            temp = read.get_thymidines()
            if len(temp) == 0:
                continue
            temp.loc[:, "sample_name"] = self.name
            temp.loc[:, "read_id"] = read._read_id + read._chrom + str(read._ref_start) + str(read._ref_end) + str(read._strand)
            temp.loc[:, "read_length"] = read.get_length()
            output.append(temp)
        return pd.concat(output)
    def __repr__(self):
        return f"<Sample name={self.name}| no. reads={len(self._reads)}"

class BrdURead():
    def __init__(self, id, chrom, ref_start, ref_end, strand):
        self._read_id = id
        self._chrom = chrom
        self._ref_start = ref_start
        self._ref_end = ref_end
        self._strand = strand
        self._thymidines = []
    def add_thymidine(self, position_on_ref, probability):
        self._thymidines.append({"chrom": self._chrom,"pos": position_on_ref, "prob_brdu": probability})
    def get_length(self):
        return abs(self._ref_end - self._ref_start)
    def get_thymidines(self):
        return pd.DataFrame(self._thymidines)
    def __repr__(self):
        return f"<BrdURead id:{self._read_id}>"

    
def read_detect_file(detect_path, limit):
    f = open(detect_path,'r')
    index = 0
    current_read = None
    for line_number, line in enumerate(f):
        #ignore the header lines
        if line[0] == '#':
                continue
        #split the line into a list by whitespace
        splitLine = line.rstrip().split()
        if line[0] == '>':
            if current_read is not None:
                    yield current_read
            index += 1
            if index > limit:
                    return
            readID = splitLine[0][1:]
            chromosome = splitLine[1]
            refStart = int(splitLine[2])
            refEnd = int(splitLine[3])
            strand = splitLine[4]
            current_read = BrdURead(
                    readID,
                    chromosome,
                    refStart,
                    refEnd,
                    strand
            )
        else:
            posOnRef = int(splitLine[0])
            probBrdU = float(splitLine[1])
            sixMerOnRef = splitLine[2]
            current_read.add_thymidine(posOnRef, probBrdU)
            #add these values to a container or do some processing here

    f.close()

def load_samples(sample_mapping):
    samples = []
    for name, paths in sample_mapping.items():
        for path in paths:
            # name with batch
            matched = re.findall(BATCH_ID, path)
            if len(matched) == 0:
                continue
            batch = matched[0]
            current_sample = Sample(f"{name}:{batch}")
            for read in read_detect_file(path, limit=100):
                current_sample.add_read(read)
            samples.append(current_sample)
    return samples

def get_thymidines(samples):
    thymidines = []
    for sample in samples:
        thymidines.append(sample.get_thymidines())
    thy_frame = pd.concat(thymidines)
    thy_frame.sample_name = thy_frame.sample_name.str.split(":").str[0]
    return thy_frame

def plot_brdu_distribution(thy_frame, ax):
    for sample in set(thy_frame.sample_name):
        ax.hist(thy_frame.query("sample_name == @sample").prob_brdu, bins=np.arange(0, 1, 0.01), label=f"{sample}", histtype="step", density=True, lw=4)
    plt.legend()
    sbn.despine()
    ax.set(title="BrdU probability distribution", xlabel="BrdU probability", ylabel="Density")

def plot_brdu_per_thymidine_distribution(thy_frame, ax):
    # get data
    thy_frame_filtered = thy_frame.query("read_length > 100")
    thy_frame_filtered.loc[:, "is_brdu"] = thy_frame_filtered.prob_brdu > 0.8
    brdu_per_read = thy_frame_filtered.groupby(["sample_name", "read_id"])\
                                      .agg({"is_brdu": [np.sum, len], "read_length": np.mean, "pos": [np.min, np.max]}).reset_index()
    brdu_per_read.loc[:, "brdu_per_thy"] = brdu_per_read.is_brdu["sum"]/brdu_per_read.is_brdu["len"]
    brdu_per_read.loc[:, "brdu_per_thy_log"] = np.log2(brdu_per_read.brdu_per_thy + 0.001)
    brdu_per_read.loc[:, "labelled"] = brdu_per_read.brdu_per_thy > 0.1
    # plot
    for sample in set(brdu_per_read.sample_name):
        ax.hist(brdu_per_read.loc[brdu_per_read.sample_name == sample, "brdu_per_thy_log"], histtype="step", density=True, lw=3, bins=50, label=sample, range=(-11, 0))
    ax.set(title="Distribution of BrdU-percentage per read [log2]", xlabel="Brdu/Thymidine [log2]", ylabel="density")
    plt.legend()

def percentage_labelled_reads(label_library_path):
    with open(label_library_path, 'rb') as handle:
        label_library = pickle.load(handle)
    return np.mean(np.array(list(label_library.values())))

def plot_fraction_labelled(label_frame, ax):
    sbn.barplot(x="sample", y="Label fraction", data=label_frame, ax=ax, color="tab:blue")
    ax.set(title="Fraction labelled")
    plt.sca(ax)
    plt.xticks(*plt.xticks(), rotation=90)
    sbn.despine()

def count_contacts(pairs_path):
    header, pairs_body = pairtools._headerops.get_header(
        pairtools._fileio.auto_open(pairs_path, "r")
    )
    count = 0 
    for i in pairs_body:
        count += 1
    return count

def plot_contact_distribution(count_frame_molten, ax):
    sbn.barplot(x="sample", y="read_number", hue="contact_type", data=count_frame_molten, ax=ax, hue_order=["all", "usable_all_reads", "cis_all_reads", "trans_all_reads", "usable_labelled", "cis_labelled", "trans_labelled", "usable_mq", "cis_mq", "trans_mq"])
    ax.set(title="Contact distribution")
    plt.sca(ax)
    plt.xticks(*plt.xticks(), rotation=90)
    sbn.despine()

# prepare data


## detect files

sample_mapping = {}

for dpath in args.detect_paths:
    run_id = Path(dpath).name.split("_")[1]
    if run_id not in sample_mapping:
        sample_mapping[run_id] = [dpath]
    else:
        sample_mapping[run_id] += [dpath]


samples = load_samples(sample_mapping)
thy_frame = get_thymidines(samples)


## label libraries

label_perc = {
    Path(path).name.split("_")[1]: percentage_labelled_reads(path) for path in args.label_library_paths
}
label_frame = pd.DataFrame(label_perc, index=["Label fraction"]).T.reset_index().rename(columns={"index": "sample"})

## contact distribution

pairs_mapping = {}

pairs_arguments = {
    "all": args.all_pairs_paths,
    "usable_all_reads": args.all_reads_cis_trans,
    "cis_all_reads": args.all_reads_cis,
    "trans_all_reads": args.all_reads_trans,
    "usable_labelled": args.labelled_cis_trans,
    "cis_labelled": args.labelled_cis,
    "trans_labelled": args.labelled_trans,
    "usable_no_c": args.no_c_cis_trans,
    "cis_no_c": args.no_c_cis,
    "trans_no_c": args.no_c_trans,
    "usable_mq": args.mq_cis_trans,
    "mq_cis": args.mq_cis,
    "mq_trans": args.mq_trans
}

for name, pairs_list in pairs_arguments.items():
    for p in pairs_list:
        run_id = Path(p).name.split("_")[1]
        if run_id in pairs_mapping:
            pairs_mapping[run_id][name] = p
        else:
            pairs_mapping[run_id] = {
                name: p
            }

# count reads

read_count = {
    sample: {
        read_type:
            count_contacts(path)
                for read_type, path in pairs_mapping[sample].items()
    }
    for sample in pairs_mapping.keys()
}

count_frame = pd.DataFrame(read_count).T.reset_index().rename(columns={"index": "sample"})
count_frame_molten = pd.melt(count_frame, id_vars="sample").rename(columns={"variable": "contact_type", "value": "read_number"})


# make plot

fig = plt.figure(constrained_layout=True)
fig.set_size_inches(14, 18)
gs = fig.add_gridspec(
    3,
    8,
    width_ratios=[1] * 8,
    height_ratios=[1, 3, 5],
    hspace=0.15,
    wspace=0.1,
)
empty = fig.add_subplot(gs[0, :])
make_title(
"Sister Pore-C report",
ax=empty,
subtitle=f"Contains information regarding distribution of labelled reads.",
)
plot_brdu_distribution(thy_frame, fig.add_subplot(gs[1, :4]))
plot_brdu_per_thymidine_distribution(thy_frame, fig.add_subplot(gs[1, 4:]))
plot_fraction_labelled(label_frame, fig.add_subplot(gs[2, :3]))
plot_contact_distribution(count_frame_molten, fig.add_subplot(gs[2, 3:]))
sbn.despine()


fig.savefig(args.output, bbox_inches="tight", pad_inches=0.5)