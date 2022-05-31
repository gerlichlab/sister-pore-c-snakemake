import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sbn
from pathlib import Path
import datetime
import argparse

# parse arguments

parser = argparse.ArgumentParser(description='Create pore-c qc plots')
parser.add_argument('--concatamer_table')
parser.add_argument('--concatamer_summary')
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


def plot_alignment_length(contact_frame, ax):
    # prepare data
    sample = contact_frame.sample(frac=0.1)
    unique_samples = sample.drop_duplicates(subset="read_name").query("contact_is_cis")
    read_lengths = unique_samples["align1_end"] - unique_samples["align1_start"]
    read_length_frame = pd.DataFrame(
        {
            "run_id": RUN_ID,
            "run_length": read_lengths
        }
    )
    # plot
    sbn.pointplot(
        x="run_id",
        y="run_length",
        data=read_length_frame,
        join=False,
        capsize=0.5,
        scale=1.5,
        color="tab:blue",
        ax=ax
    )
    sbn.despine()
    ax.set(xlabel="Run id", ylabel="Alignment length [bp]", title="Digested fragment size")

def plot_absolute_read_distribution(summary_table, ax):
    # Prepare data
    read_tables = summary_table.query("section in ['reads', 'read_length']")
    read_tables.loc[:, "value"] = read_tables.value.astype(float)
    # Plot
    sbn.barplot(x="level_2", y="value", data=read_tables.query("level_1 == 'count'"), ax=ax, color="tab:blue")
    ax.set(title="Read distribution abs", xlabel="Read-type", ylabel="Count")

def plot_read_length(summary_table, ax):
    # prepare data
    read_tables = summary_table.query("section in ['reads', 'read_length']")
    read_tables.loc[:, "value"] = read_tables.value.astype(float)
    # Plot
    sbn.barplot(x="level_2", y="value", data=read_tables.query("level_1 == 'N50'"), ax=ax, color="tab:blue")
    ax.set(title="Read length", xlabel="Read-type", ylabel="Read length")

def plot_contact_density(summary_table, ax):
    # prepare data
    density_table = summary_table.query("section in ['density']")
    density_table.loc[:, "value"] = density_table.value.astype("float")
    # plot
    sbn.barplot(x="run_id", y="value", data=density_table.query("level_2 == 'all'"), ax=ax, color="tab:blue")
    ax.set(title="Density of contacts/Gbp", xlabel="", ylabel="Contacts/Gbp")

def plot_contact_distribution(summary_table, ax):
    # prepare data
    concat_table = summary_table.query("section in ['concatemer_order']")
    concat_table.loc[:, "value"] = concat_table.value.astype(float)
    # plot
    sbn.barplot(x="level_2", y="value", data=concat_table.query("level_1 == 'perc'"), order=["2","3", "4", "5-6", "7-11", "12-21", "22-50", "gt_50"], ax=ax, color="tab:blue")
    ax.set(title="Concatamer distribution", ylabel="percentage", xlabel="Concatamer order")

def plot_contact_type(summary_table, ax):
    # prepare data
    contact_table = summary_table.query("section in ['contacts']")
    contact_table.loc[:, "value"] = contact_table.value.astype("float")
    # plot
    sbn.barplot(x="level_2", y="value", data=contact_table.query("level_0 == 'total' and level_1 == 'perc'"), ax=ax, color="tab:blue")
    ax.set(xlabel="Read-type", ylabel="Total contacts [%]")
    plt.sca(ax)
    plt.xticks(*plt.xticks(), rotation=90)

def plot_qc_report(contact_frame, summary_table):
    fig = plt.figure(constrained_layout=True)
    fig.set_size_inches(14, 18)
    gs = fig.add_gridspec(
        5,
        8,
        width_ratios=[1] * 8,
        height_ratios=[1, 3, 3, 3, 3],
        hspace=0.15,
        wspace=0.1,
    )
    empty = fig.add_subplot(gs[0, :])
    make_title(
        "Pore-C QC report",
        ax=empty,
        subtitle=f"Contains information regarding fragment length, contact density and concatamer order distribution.",
    )
    plot_alignment_length(contact_frame, fig.add_subplot(gs[1, :4]))
    plot_absolute_read_distribution(summary_table, fig.add_subplot(gs[1, 4:]))
    plot_read_length(summary_table, fig.add_subplot(gs[2, :4]))
    plot_contact_type(summary_table, fig.add_subplot(gs[2, 4:]))
    plot_contact_density(summary_table, fig.add_subplot(gs[3, :4]))
    plot_contact_distribution(summary_table, fig.add_subplot(gs[3, 4:]))
    sbn.despine()
    return fig

# Constants

DATE = r"(\d){4}-(\d){2}-(\d){2}"
RUN_NAME = r"(run(\d){2})"
RUN_ID = Path(args.concatamer_table).name.split("_")[1]
plt.rcParams["font.size"] = 15


# Distribution of read-lengths


contact_frame = pd.read_parquet(args.concatamer_table)
summary_table = pd.read_csv(args.concatamer_summary)
summary_table.loc[:, "run_id"] = RUN_ID


f = plot_qc_report(contact_frame, summary_table)
f.savefig(args.output, bbox_inches="tight", pad_inches=0.5)

exit()

# Load reports for different stages in the pipeline



## Contact information
contact_table = summary_table.query("section in ['contacts']")
contact_table.loc[:, "value"] = contact_table.value.astype("float")
contact_table_grouped = contact_table.groupby(["run_id", "section", "level_0", "level_1", "level_2"]).sum().reset_index()


f, ax = plt.subplots(1, 3)
plt.subplots_adjust(wspace=0.3)
sbn.barplot(x="level_2", y="value", data=contact_table_grouped.query("level_0 == 'direct' and level_1 == 'count'"), ax=ax[0], hue="run_id")
ax[0].set(xlabel="Read-type", ylabel="Direct Contacts")
sbn.barplot(x="level_2", y="value", data=contact_table_grouped.query("level_0 == 'indirect' and level_1 == 'count'"), ax=ax[1], hue="run_id")
ax[1].set(xlabel="Read-type", ylabel="Count")
sbn.barplot(x="level_2", y="value", data=contact_table_grouped.query("level_0 == 'total' and level_1 == 'count'"), ax=ax[2], hue="run_id")
ax[2].set(xlabel="Read-type", ylabel="Count")
f.set_size_inches(30, 5)
for axis in ax:
    plt.sca(axis)
    plt.xticks(*plt.xticks(), rotation=90)
for i in range(3):
    ax[i].legend(
            loc="upper center",
            bbox_to_anchor=(0.5, 1.3),
            ncol=3,
            fancybox=True,
            shadow=True,
    )
sbn.despine()
plt.show()


f, ax = plt.subplots(1, 3)
plt.subplots_adjust(wspace=0.3)
sbn.barplot(x="level_2", y="value", data=contact_table.query("level_0 == 'direct' and level_1 == 'perc'"), ax=ax[0], hue="run_id")
ax[0].set(xlabel="Read-type", ylabel="Direct contacts [%]")
sbn.barplot(x="level_2", y="value", data=contact_table.query("level_0 == 'indirect' and level_1 == 'perc'"), ax=ax[1], hue="run_id")
ax[1].set(xlabel="Read-type", ylabel="Indirect contacts [%]")
sbn.barplot(x="level_2", y="value", data=contact_table.query("level_0 == 'total' and level_1 == 'perc'"), ax=ax[2], hue="run_id")
ax[2].set(xlabel="Read-type", ylabel="Total contacts [%]")
f.set_size_inches(30, 7)
for axis in ax:
    plt.sca(axis)
    plt.xticks(*plt.xticks(), rotation=90)
for i in range(3):
    ax[i].legend(
            loc="upper center",
            bbox_to_anchor=(0.5, 1.2),
            ncol=3,
            fancybox=True,
            shadow=True,
    )
sbn.despine()
plt.show()


## Contact density
density_table = summary_table.query("section in ['density']")
density_table.loc[:, "value"] = density_table.value.astype("float")


f, ax = plt.subplots()
sbn.barplot(x="run_id", y="value", data=density_table.query("level_2 == 'all'"), ax=ax, order=density_table.query("level_2 == 'all'").sort_values("value", ascending=False).run_id)
ax.set(title="Density of contacts/Gbp", xlabel="", ylabel="Contacts/Gbp")
f.set_size_inches(13, 5)
sbn.despine()
plt.show()

# Concatamer order

concat_table = summary_table.query("section in ['concatemer_order']")
concat_table.loc[:, "value"] = concat_table.value.astype(float)
concat_table.query("level_1 == 'count'").sum()

## Add singletons

pore_c_frame.loc[:, "1"] = (pore_c_frame.singleton / 50000) * 100 # batch size was 50k
sing_frame = pore_c_frame[["run_id", "1"]].groupby("run_id").mean().reset_index().melt(id_vars="run_id").rename(columns={"variable": "level_2"})
sing_frame.loc[:, "level_1"] = "perc"
concat_w_singletons = pd.concat((concat_table, sing_frame))


f, ax = plt.subplots()
sbn.barplot(x="level_2", y="value", data=concat_w_singletons.query("level_1 == 'perc'"), order=["1", "2","3", "4", "5-6", "7-11", "12-21", "22-50", "gt_50"], ax=ax, hue="run_id")
ax.set(title="Concatamer distribution", ylabel="percentage", xlabel="Concatamer order")
f.set_size_inches(14, 5)
sbn.despine()
plt.show()