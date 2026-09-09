import argparse
import pandas as pd
import matplotlib.pyplot as plt
import met_brewer


FONT_SIZE = 7
FIG_WIDTH_IN = 160 / 25.4
FIG_HEIGHT_IN = 30 / 25.4


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--admixture-afts-pkl",
        required=True,
        help=("Path to admixture_afts.pkl."),
    )
    parser.add_argument(
        "--output-pdf",
        required=True,
        help=("Path to output PDF file."),
    )
    args = parser.parse_args()

    simulations_df = pd.read_pickle(args.admixture_afts_pkl)

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

    props = [0.1, 0.2, 0.5]
    titles = {
        0.1: "Admixture proportion = 0.1",
        0.2: "Admixture proportion = 0.2",
        0.5: "Admixture proportion = 0.5",
    }

    for i, (ax, prop) in enumerate(zip(axes, props)):
        subset = simulations_df[
            (simulations_df["S1"] == 0.0)
            & (simulations_df["S2"] == 1.0)
        ]
        subset = subset[subset["PROP"] == prop]
        subset["ticks"] = subset["ticks"] - 150000
        subset = subset.sort_values("seed").head(10)

        colors = met_brewer.met_brew(name="Redon", n=10)

        for color, (_, row) in zip(colors, subset.iterrows()):
            ticks = row["ticks"]
            host_A_values = row["host_A_series"]
            symb_A_values = row["symb_A_series"]
            seed = row["seed"]

            ax.plot(
                ticks,
                host_A_values,
                color=color,
                linestyle="-",
                alpha=0.8
            )

            ax.plot(
                ticks,
                symb_A_values,
                color=color,
                linestyle="--",
                alpha=0.8
            )
        ax.plot([0,1], [0, 0],
        color="grey",
        linestyle="--",
        alpha=0.8,
        label=f"Symbiont trajectory (a)",)

        ax.plot([0,1], [0, 0],
                color="grey",
                linestyle="-",
                alpha=0.8,
                label=f"Host trajectory (A)",)
        ax.set_xlim(1, 100)
        ax.set_ylim(0, 1)
        ax.spines["top"].set_alpha(0)
        ax.spines["right"].set_alpha(0)
        ax.set_title(titles[prop], fontsize=FONT_SIZE)
        ax.tick_params(axis="both", labelsize=FONT_SIZE)

        if i == 0:
            ax.set_xlabel("Generation", fontsize=FONT_SIZE)
            ax.set_ylabel("Allele frequency", fontsize=FONT_SIZE)
            ax.legend(loc="upper right")
        else:
            ax.set_xlabel("", fontsize=FONT_SIZE)
            ax.set_ylabel("", fontsize=FONT_SIZE)

    fig.savefig(args.output_pdf)
    plt.close(fig)


if __name__ == "__main__":
    main()