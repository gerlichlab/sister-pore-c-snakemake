import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sbn
from pathlib import Path
import datetime
import argparse

# parse arguments

parser = argparse.ArgumentParser(description="Create pore-c qc plots")
parser.add_argument("--concatamer_tables", nargs="+")
parser.add_argument("--concatamer_summaries", nargs="+")
parser.add_argument("--output")
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
    unique_samples = contact_frame.drop_duplicates(subset="read_name").query(
        "contact_is_cis"
    )
    read_lengths = unique_samples["align1_end"] - unique_samples["align1_start"]
    read_length_frame = pd.DataFrame(
        {"run_id": unique_samples.run_id, "run_length": read_lengths}
    )
    # plot
    sbn.pointplot(
        x="run_id",
        y="run_length",
        data=read_length_frame,
        join=False,
        capsize=0.5,
        scale=1.5,
        ax=ax,
    )
    sbn.despine()
    ax.set(
        xlabel="Run id", ylabel="Alignment length [bp]", title="Digested fragment size"
    )


def plot_absolute_read_distribution(summary_table, ax):
    # Prepare data
    read_tables = summary_table.query("section in ['reads', 'read_length']")
    read_tables.loc[:, "value"] = read_tables.value.astype(float)
    # Plot
    sbn.barplot(
        x="level_2",
        y="value",
        hue="run_id",
        data=read_tables.query("level_1 == 'count'"),
        ax=ax,
    )
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.2),
        ncol=3,
        fancybox=True,
        shadow=True,
    )
    ax.set(title="Read distribution abs", xlabel="Read-type", ylabel="Count")


def plot_read_length(summary_table, ax):
    # prepare data
    read_tables = summary_table.query("section in ['reads', 'read_length']")
    read_tables.loc[:, "value"] = read_tables.value.astype(float)
    # Plot
    sbn.barplot(
        x="level_2",
        y="value",
        hue="run_id",
        data=read_tables.query("level_1 == 'N50'"),
        ax=ax,
    )
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.2),
        ncol=3,
        fancybox=True,
        shadow=True,
    )
    ax.set(title="Read length", xlabel="Read-type", ylabel="Read length")


def plot_contact_density(summary_table, ax):
    # prepare data
    density_table = summary_table.query("section in ['density']")
    density_table.loc[:, "value"] = density_table.value.astype("float")
    # plot
    sbn.barplot(
        x="run_id", y="value", data=density_table.query("level_2 == 'all'"), ax=ax
    )
    ax.set(title="Density of contacts/Gbp", xlabel="", ylabel="Contacts/Gbp")


def plot_contact_distribution(summary_table, ax):
    # prepare data
    concat_table = summary_table.query("section in ['concatemer_order']")
    concat_table.loc[:, "value"] = concat_table.value.astype(float)
    # plot
    sbn.barplot(
        x="level_2",
        y="value",
        hue="run_id",
        data=concat_table.query("level_1 == 'perc'"),
        order=["2", "3", "4", "5-6", "7-11", "12-21", "22-50", "gt_50"],
        ax=ax,
    )
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.2),
        ncol=3,
        fancybox=True,
        shadow=True,
    )
    ax.set(
        title="Concatamer distribution", ylabel="percentage", xlabel="Concatamer order"
    )


def plot_contact_type(summary_table, ax):
    # prepare data
    contact_table = summary_table.query("section in ['contacts']")
    contact_table.loc[:, "value"] = contact_table.value.astype("float")
    # plot
    sbn.barplot(
        x="level_2",
        y="value",
        hue="run_id",
        data=contact_table.query("level_0 == 'total' and level_1 == 'perc'"),
        ax=ax,
    )
    ax.set(xlabel="Read-type", ylabel="Total contacts [%]")
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.2),
        ncol=3,
        fancybox=True,
        shadow=True,
    )
    plt.sca(ax)
    plt.xticks(*plt.xticks(), rotation=90)


def plot_qc_report(contact_frame, summary_table, page=None):
    fig = plt.figure(constrained_layout=True)
    fig.set_size_inches(15, 20)
    gs = fig.add_gridspec(
        5,
        8,
        width_ratios=[1] * 8,
        height_ratios=[1, 3, 3, 3, 3],
        hspace=0.25,
        wspace=0.15,
    )
    empty = fig.add_subplot(gs[0, :])
    if page is None:
        make_title(
            "Pore-C QC report",
            ax=empty,
            subtitle=f"Contains information regarding fragment length, contact density and concatamer order distribution.",
        )
    else:
        make_title(
            f"Pore-C QC report | page {page}",
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


def get_run_id(path):
    return Path(path).name.split("_")[1]


# Set constants and plotting parameters

DATE = r"(\d){4}-(\d){2}-(\d){2}"
plt.rcParams["font.size"] = 12


# read input files and associate run id

# contacts
contact_frames = []
for contact_frame_path in args.concatamer_tables:
    run_id = get_run_id(contact_frame_path)
    temp_contact_frame = pd.read_parquet(contact_frame_path)
    temp_contact_frame.loc[:, "run_id"] = run_id
    contact_frames.append(temp_contact_frame.sample(frac=0.1))

contact_frame = pd.concat(contact_frames)

# summary tables

summary_tables = []
for summary_path in args.concatamer_summaries:
    run_id = get_run_id(summary_path)
    temp_summary_frame = pd.read_csv(summary_path)
    temp_summary_frame.loc[:, "run_id"] = run_id
    summary_tables.append(temp_summary_frame)

summary_table = pd.concat(summary_tables)

# if number of samples is below 5, plot single page, if not, plot multiple pages

samples = sorted(list(set(summary_table.run_id)))
number_samples = len(samples)

if number_samples < 5:
    f = plot_qc_report(contact_frame, summary_table)
    f.savefig(args.output, bbox_inches="tight", pad_inches=0.5)
else:
    with PdfPages(args.output) as pdf:
        number_pages = number_samples // 5 + 1
        for i in range(0, number_pages):
            current_samples = samples[i * 5 : (i + 1) * 5]
            subset_contact_frame = contact_frame.query("run_id in @current_samples")
            subset_summary_table = summary_table.query("run_id in @current_samples")
            f = plot_qc_report(subset_contact_frame, subset_summary_table, page=i)
            pdf.savefig(f, bbox_inches="tight", pad_inches=0.5)
