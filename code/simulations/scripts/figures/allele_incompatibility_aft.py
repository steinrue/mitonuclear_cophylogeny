import argparse
import pandas as pd
import matplotlib.pyplot as plt
import met_brewer


FONT_SIZE = 7
FIG_WIDTH_IN = 80 / 25.4
FIG_HEIGHT_IN = 30 / 25.4


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--allele-incompatability-afts-pkl",
        required=True,
        help=("Path to allele_incompatability_afts.pkl."),
    )
    parser.add_argument(
        "--output-pdf",
        required=True,
        help=("Path to output PDF file."),
    )
    args = parser.parse_args()

    allele_incompatability_freqs_df = pd.read_pickle(
        args.allele_incompatability_afts_pkl
    )

    subset = allele_incompatability_freqs_df[
        (allele_incompatability_freqs_df["S1"] == 0.0)
        & (allele_incompatability_freqs_df["S2"] == 1)
    ]
    subset = subset[
        subset["symb_A_series"].apply(
            lambda arr: (arr[0] > 0.375) and (arr[0] < 0.625)
        )
    ]
    subset = subset.sort_values("seed")
    passing_seeds = list(subset["seed"])
    subset["ticks"] = subset["ticks"] - 150009

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

    fig, ax = plt.subplots(
        figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
        layout="constrained",
        dpi=600,
    )

    colors = met_brewer.met_brew(name="Redon", n=len(passing_seeds))

    for color, (_, row) in zip(colors, subset.iterrows()):
        ticks = row["ticks"]
        host_A_values = row["host_A_series"]
        symb_A_values = row["symb_A_series"]

        ax.plot(
            ticks,
            host_A_values,
            color=color,
            linestyle="-",
            alpha=0.8,
        )

        ax.plot(
            ticks,
            symb_A_values,
            color=color,
            linestyle="--",
            alpha=0.8,
        )

    ax.plot(
        [0, 1],
        [0, 0],
        color="grey",
        linestyle="--",
        alpha=0.8,
        label="Symbiont trajectory (a)",
    )

    ax.plot(
        [0, 1],
        [0, 0],
        color="grey",
        linestyle="-",
        alpha=0.8,
        label="Host trajectory (A)",
    )

    ax.set_xlim(0, 100)
    ax.set_xlabel("Generation", fontsize=FONT_SIZE)
    ax.set_ylabel("Allele frequency", fontsize=FONT_SIZE)
    ax.spines["top"].set_alpha(0)
    ax.spines["right"].set_alpha(0)
    ax.set_ylim(0, 1)
    ax.legend(loc="center right",fontsize=FONT_SIZE,)

    fig.savefig(args.output_pdf)
    plt.close(fig)


if __name__ == "__main__":
    main()