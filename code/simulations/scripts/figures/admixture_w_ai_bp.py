import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from met_brewer import met_brew

FONT_SIZE = 7
FIG_WIDTH_IN = 100 / 25.4
FIG_HEIGHT_IN = 40 / 25.4

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--admixture-pkl",
        required=True,
        help=("Path to admixture.pkl."),
    )
    parser.add_argument(
        "--output-pdf",
        required=True,
        help=("Path to output PDF file."),
    )
    args = parser.parse_args()

    df1 = pd.read_pickle(args.admixture_pkl)
    df1["global_mean_zscore"] = -df1["global_mean_zscore"]
    df1 = df1[(df1["ST"] != 0)]
    df1 = df1[(df1["S1"] == 0) & (df1["S2"] == 1)]

    plot_df = df1.copy()
    prop_order = sorted(plot_df["PROP"].unique())
    plot_df["PROP"] = pd.Categorical(plot_df["PROP"], categories=prop_order, ordered=True)

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

    fig, axs = plt.subplots(
        figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
        layout="constrained",
        dpi=600,
    )

    grouped = sorted(plot_df["ST"].unique())
    colors = met_brew(
        name="Hokusai2",
        n=len(grouped),
        brew_type="continuous",
    )
    palette = dict(zip(grouped, colors))

    sns.boxplot(
        data=plot_df,
        x="PROP",
        y="global_mean_zscore",
        hue="ST",
        palette=palette,
        ax=axs,
        linewidth=0.75,
        showfliers=True,
        flierprops=dict(
            marker=".",
            markersize=1,
            markerfacecolor="black",
            markeredgecolor="black",
        ),
    )

    axs.legend(
        loc="upper right",
        bbox_to_anchor=(1.2, 1),
        borderaxespad=0,
        title="Generations\nafter admixture",
        fontsize=FONT_SIZE,
        title_fontsize=FONT_SIZE,
        frameon=True,
        edgecolor="none",
        framealpha=0,
    )

    axs.axhline(0, color="black", linestyle="--", linewidth=1)

    axs.tick_params(axis="both", which="both", length=0, labelsize=FONT_SIZE)
    axs.set_xlabel(r"Admixture proportion", fontsize=FONT_SIZE)
    axs.set_ylabel(r"Cophylogeny score", fontsize=FONT_SIZE)
    axs.spines["top"].set_visible(False)
    axs.spines["right"].set_visible(False)

    fig.savefig(args.output_pdf)
    plt.close(fig)


if __name__ == "__main__":
    main()