#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

FONT_SIZE = 7
FIG_WIDTH_IN = 120 / 25.4
FIG_HEIGHT_IN = 40 / 25.4
MAIN_COLOR = "darkgreen"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tips50", required=True, help="Path to tips_50_perms_100000.pkl")
    parser.add_argument("--tips100", required=True, help="Path to tips_100_perms_100000.pkl")
    parser.add_argument("--tips200", required=True, help="Path to tips_200_perms_100000.pkl")
    parser.add_argument("--output", required=True, help="Path to output PDF/PNG file")
    args = parser.parse_args()

    files = {
        "50 haplotypes": args.tips50,
        "100 haplotypes": args.tips100,
        "200 haplotypes": args.tips200,
    }

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
        3,
        figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
        layout="constrained",
        dpi=600,
    )

    for ax, (title, path) in zip(axes, files.items()):
        df = pd.read_pickle(path)

        rep_cols = [c for c in df.columns if c not in ["expected_p", "lower", "upper"]]
        if len(rep_cols) == 0:
            raise ValueError(f"No replicate columns found in {path}. Expected columns other than 'expected_p', 'lower', 'upper'.")

        p_mat = df[rep_cols].to_numpy()
        p_mat_log = -np.log10(p_mat)

        medians = np.nanmedian(p_mat_log, axis=1)
        stds = np.nanstd(p_mat_log, axis=1, ddof=1)

        exp = -np.log10(df["expected_p"].to_numpy())
        # lower = df["lower"].to_numpy()
        # upper = df["upper"].to_numpy()
        # ax.fill_between(exp, lower, upper, color="k", alpha=0.1, rasterized=True)

        ax.plot([0, 6], [0, 6], linestyle="--", color="grey", linewidth=1)

        ax.scatter(
            exp,
            medians,
            alpha=0.8,
            s=8,
            color=MAIN_COLOR,
            edgecolor="none",
            rasterized=True,
        )

        ax.fill_between(
            exp,
            medians - stds,
            medians + stds,
            color=MAIN_COLOR,
            alpha=0.2,
            rasterized=True,
        )

        ax.set_xlim(0, 6)
        ax.set_ylim(0, 6)
        ax.spines["top"].set_alpha(0)
        ax.spines["right"].set_alpha(0)
        ax.set_title(title, fontsize=FONT_SIZE)

        ax.tick_params(axis="both", labelsize=FONT_SIZE, length=2)
        ax.grid(False)

        axins = inset_axes(
            ax,
            width="30%",
            height="30%",
            loc="lower right",
        )

        raw_exp = np.sort(10 ** (-exp))
        raw_medians = np.sort(10 ** (-medians))

        axins.plot([0, 1], [0, 1], linestyle="--", color="grey", linewidth=1)

        axins.scatter(
            raw_exp[:len(raw_medians)],
            raw_medians,
            s=3,
            alpha=0.8,
            color=MAIN_COLOR,
            edgecolor="none",
            rasterized=True,
        )

        axins.set_xlim(0, 1)
        axins.set_ylim(0, 1)
        axins.tick_params(axis="y", which="both", left=False)
        axins.yaxis.tick_left()
        axins.yaxis.set_label_position("left")
        axins.set_xticks([])
        axins.tick_params(axis="both", labelsize=FONT_SIZE, length=2)
        axins.grid(False)

    for ax in axes:
        ax.set_xlabel("")

    axes[0].set_ylabel("Permutation estimation\n" r"$-\log_{10}(p)$", fontsize=FONT_SIZE)
    axes[1].set_ylabel("")
    axes[2].set_ylabel("")

    fig.supxlabel(r"Expected $-\log_{10}(p)$", fontsize=FONT_SIZE)

    fig.savefig(args.output, dpi=600)
    plt.close(fig)

if __name__ == "__main__":
    main()