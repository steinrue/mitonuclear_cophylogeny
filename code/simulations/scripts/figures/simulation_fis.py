import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import met_brewer
from met_brewer import met_brew
from scipy.stats import norm

FONT_SIZE = 7
FIG_WIDTH_IN = 160 / 25.4
FIG_HEIGHT_IN = 40 / 25.4

def plot_unstructured(ax, df):
    # Prepare data for boxplot
    groups = sorted(df["K"].unique())
    data = [df[df["K"] == g]["fis"] for g in groups]

    colors = met_brew(name="Homer2", brew_type="continuous", n=6)[3:]

    bplot = ax.boxplot(
        data,
        patch_artist=True,
        tick_labels=groups,
        showfliers=False,
        medianprops={'color': 'black'}
    )

    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

    ax.set_ylabel(r"$F_{IS}$", fontsize=FONT_SIZE)
    ax.set_xlabel("Population size", fontsize=FONT_SIZE)
    ax.axhline(0.0, ls="--", color="grey", alpha=0.8, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("Unstructured population", fontsize=FONT_SIZE)
    ax.tick_params(axis="y", labelsize=FONT_SIZE)
    ax.tick_params(axis="x", labelsize=FONT_SIZE)


def plot_migration(ax, df):
    groups = sorted(df["M"].unique())
    data = [df[df["M"] == g]["fis"] for g in groups]

    ordered_colors = met_brewer.met_brew(
        name="OKeeffe1",
        n=(len(groups) + 2) * 2,
        brew_type="continuous",
    )
    colors = ordered_colors[:len(groups)]

    bplot = ax.boxplot(
        data,
        patch_artist=True,
        tick_labels=groups,
        showfliers=False,
        medianprops={'color': 'black'}
    )

    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

    ax.set_ylabel(r"", fontsize=FONT_SIZE)
    ax.set_xlabel("Migration rate", fontsize=FONT_SIZE)
    ax.axhline(0.0, ls="--", color="grey", alpha=0.8, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("Migration", fontsize=FONT_SIZE)
    ax.tick_params(axis="y", labelsize=FONT_SIZE)
    ax.tick_params(axis="x", labelsize=FONT_SIZE)


def plot_pop_split(ax, df):
    groups = sorted(df["ST"].unique())
    data = [df[df["ST"] == g]["fis"] for g in groups]

    ordered_colors = met_brewer.met_brew(
        name="Benedictus",
        n=(len(groups) + 2) * 2,
        brew_type="continuous",
    )
    colors = ordered_colors[:len(groups)]

    bplot = ax.boxplot(
        data,
        patch_artist=True,
        tick_labels=groups,
        showfliers=False,
        medianprops={'color': 'black'}
    )

    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

    ax.set_ylabel(r"", fontsize=FONT_SIZE)
    ax.set_xlabel("Generations after split", fontsize=FONT_SIZE)
    ax.axhline(0.0, ls="--", color="grey", alpha=0.8, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("Population split", fontsize=FONT_SIZE)
    ax.tick_params(axis="y", labelsize=FONT_SIZE)
    ax.tick_params(axis="x", labelsize=FONT_SIZE)


def plot_admixture(ax, df):
    df = df[df["PROP"] == 0.5]
    df = df[df["ST"] != 0]
    df = df[df["ST"] != 100]

    groups = sorted(df["ST"].unique())
    data = [df[df["ST"] == g]["fis"] for g in groups]

    colors = met_brewer.met_brew(
        name="Hokusai2",
        n=len(groups),
        brew_type="continuous",
    )

    bplot = ax.boxplot(
        data,
        patch_artist=True,
        tick_labels=groups,
        showfliers=False,
        medianprops={'color': 'black'}
    )

    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

    ax.set_ylabel(r"", fontsize=FONT_SIZE)
    ax.set_xlabel("Generations after admixture", fontsize=FONT_SIZE)
    ax.axhline(0.0, ls="--", color="grey", alpha=0.8, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("Admixture", fontsize=FONT_SIZE)
    ax.tick_params(axis="y", labelsize=FONT_SIZE)
    ax.tick_params(axis="x", labelsize=FONT_SIZE)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fis-tm-pkl",
        required=True,
        help="Path to fis_tm.pkl.",
    )
    parser.add_argument(
        "--fis-mig-pkl",
        required=True,
        help="Path to fis_mig.pkl.",
    )
    parser.add_argument(
        "--fis-pop-split-pkl",
        required=True,
        help="Path to fis_pop_split.pkl.",
    )
    parser.add_argument(
        "--fis-admixture-pkl",
        required=True,
        help="Path to fis_admixture.pkl.",
    )
    parser.add_argument(
        "--output-pdf",
        required=True,
        help="Path to output PDF file.",
    )
    args = parser.parse_args()

    df_tm = pd.read_pickle(args.fis_tm_pkl)
    df_mig = pd.read_pickle(args.fis_mig_pkl)
    df_pop_split = pd.read_pickle(args.fis_pop_split_pkl)
    df_admixture = pd.read_pickle(args.fis_admixture_pkl)

    plt.rcParams.update({
    'font.size': FONT_SIZE,
    'text.usetex': False,
    'font.family': 'sans-serif',
    'font.sans-serif': 'Liberation Sans',
    'mathtext.fontset': 'stixsans',
    'axes.unicode_minus': False,
    'axes.formatter.use_mathtext': True,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none'
    })

    fig, axes = plt.subplots(
        1,
        4,
        figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
        layout="constrained",
        sharex=False,
        sharey=False,
        dpi=600,
    )

    plot_unstructured(axes[0], df_tm)
    plot_migration(axes[1], df_mig)
    plot_pop_split(axes[2], df_pop_split)
    plot_admixture(axes[3], df_admixture)

    # Removed fig.supxlabel(Fis) since it is now on the Y-axes of individual plots
    
    plt.savefig(args.output_pdf, format="pdf", bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()