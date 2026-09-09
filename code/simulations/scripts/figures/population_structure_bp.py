import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from met_brewer import met_brew

FONT_SIZE = 7
FIG_WIDTH_IN = 170 / 25.4
FIG_HEIGHT_IN = 50 / 25.4

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--migration-pkl", required=True)
    parser.add_argument("--symp-migration-pkl", required=True)
    parser.add_argument("--population-split-pkl", required=True)
    parser.add_argument("--symp-population-split-pkl", required=True)
    parser.add_argument("--admixture-pkl", required=True)
    parser.add_argument("--output-pdf", required=True)
    args = parser.parse_args()

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
        1, 3,
        figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
        sharey=True,
        layout="constrained",
        dpi=600,
    )

    df1 = pd.read_pickle(args.migration_pkl)
    df1["global_mean_zscore"] = -df1["global_mean_zscore"]
    df2 = pd.read_pickle(args.symp_migration_pkl)
    df2["global_mean_zscore"] = -df2["global_mean_zscore"]
    vertical = df1[df1["PM"] == 1.0].copy()
    vertical["Transmission"] = "Vertical"
    allo = df1[df1["PR"] == 1.0].copy()
    allo["Transmission"] = "Horizontal\n(allopatric)"
    symp = df2[df2["PR"] == 1.0].copy()
    symp["Transmission"] = "Horizontal\n(sympatric)"
    plot_df = pd.concat([vertical, allo, symp], axis=0)
    plot_df["Transmission"] = pd.Categorical(
        plot_df["Transmission"],
        categories=["Vertical", "Horizontal\n(allopatric)", "Horizontal\n(sympatric)"],
        ordered=True,
    )
    migration_rates = sorted(plot_df["M"].unique())
    colors = met_brew(name="OKeeffe1", n=(len(migration_rates) + 2) * 2, brew_type="continuous")[:len(migration_rates)]
    palette = dict(zip(migration_rates, colors))
    sns.boxplot(
        data=plot_df, x="Transmission", y="global_mean_zscore", hue="M",
        palette=palette, ax=axs[0], linewidth=0.75, showfliers=True,
        flierprops=dict(marker=".", markersize=1, markerfacecolor="black", markeredgecolor="black")
    )
    axs[0].legend(loc="upper right", borderaxespad=0, title="Migration rate", fontsize=FONT_SIZE, title_fontsize=FONT_SIZE, frameon=False, edgecolor="none", framealpha=0)
    axs[0].set_xlabel(r"Transmission mode", fontsize=FONT_SIZE)
    axs[0].set_ylabel(r"Cophylogeny score", fontsize=FONT_SIZE)

    df1 = pd.read_pickle(args.population_split_pkl)
    df1["global_mean_zscore"] = -df1["global_mean_zscore"]
    df2 = pd.read_pickle(args.symp_population_split_pkl)
    df2["global_mean_zscore"] = -df2["global_mean_zscore"]
    df1 = df1[df1["ST"] != 0]
    df2 = df2[df2["ST"] != 0]
    vertical = df1[df1["PM"] == 1.0].copy()
    vertical["Transmission"] = "Vertical"
    allo = df1[df1["PR"] == 1.0].copy()
    allo["Transmission"] = "Horizontal\n(allopatric)"
    symp = df2[df2["PR"] == 1.0].copy()
    symp["Transmission"] = "Horizontal\n(sympatric)"
    plot_df = pd.concat([vertical, allo, symp], axis=0)
    plot_df["Transmission"] = pd.Categorical(
        plot_df["Transmission"],
        categories=["Vertical", "Horizontal\n(allopatric)", "Horizontal\n(sympatric)"],
        ordered=True,
    )
    split_times = sorted(plot_df["ST"].unique())
    colors = met_brew(name="Benedictus", n=(len(split_times) + 2) * 2, brew_type="continuous")[:len(split_times)]
    palette = dict(zip(split_times, colors))
    sns.boxplot(
        data=plot_df, x="Transmission", y="global_mean_zscore", hue="ST",
        palette=palette, ax=axs[1], linewidth=0.75, showfliers=True,
        flierprops=dict(marker=".", markersize=1, markerfacecolor="black", markeredgecolor="black")
    )
    axs[1].legend(loc="upper right", borderaxespad=0, title="Generations\nafter split", fontsize=FONT_SIZE, title_fontsize=FONT_SIZE, frameon=False, edgecolor="none", framealpha=0)
    axs[1].set_xlabel(r"Transmission mode", fontsize=FONT_SIZE)
    axs[1].set_ylabel("")

    df1 = pd.read_pickle(args.admixture_pkl)
    df1["global_mean_zscore"] = -df1["global_mean_zscore"]
    df1 = df1[(df1["ST"] != 0) & (df1["ST"] != 100)]
    df1 = df1[(df1["S1"] == 0) & (df1["S2"] == 0)]
    plot_df = df1.copy()
    prop_order = sorted(plot_df["PROP"].unique())
    plot_df["PROP"] = pd.Categorical(plot_df["PROP"], categories=prop_order, ordered=True)
    grouped = sorted(plot_df["ST"].unique())
    colors = met_brew(name="Hokusai2", n=len(grouped), brew_type="continuous")
    palette = dict(zip(grouped, colors))
    sns.boxplot(
        data=plot_df, x="PROP", y="global_mean_zscore", hue="ST",
        palette=palette, ax=axs[2], linewidth=0.75, showfliers=True,
        flierprops=dict(marker=".", markersize=1, markerfacecolor="black", markeredgecolor="black")
    )
    axs[2].legend(loc="upper right", bbox_to_anchor=(1.15, 1), borderaxespad=0, title="Generations\nafter admixture", fontsize=FONT_SIZE, title_fontsize=FONT_SIZE, frameon=True, edgecolor="none", framealpha=0)
    axs[2].set_xlabel(r"Admixture proportion", fontsize=FONT_SIZE)
    axs[2].set_ylabel("")

    for ax in axs:
        ax.axhline(0, color="black", linestyle="--", linewidth=1)
        ax.tick_params(axis="both", which="both", length=0, labelsize=FONT_SIZE)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.savefig(args.output_pdf)
    plt.close(fig)

if __name__ == "__main__":
    main()