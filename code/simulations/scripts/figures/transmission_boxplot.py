import argparse
import matplotlib.gridspec as gridspec
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from met_brewer import met_brew


FONT_SIZE = 7
FIG_WIDTH_IN = 120 / 25.4
FIG_HEIGHT_IN = 30 / 25.4


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--transmission-w-filter-pkl",
        required=True,
        help="Path to transmission_w_filter.pkl (filtered transmission data).",
    )
    parser.add_argument(
        "--transmission-pkl",
        required=True,
        help="Path to transmission.pkl (full transmission data).",
    )
    parser.add_argument(
        "--output-pdf",
        required=True,
        help="Path to output PDF file for the figure.",
    )
    args = parser.parse_args()

    df_filtered = pd.read_pickle(args.transmission_w_filter_pkl)
    df_filtered["global_mean_zscore"] = -df_filtered["global_mean_zscore"]
    df_all = pd.read_pickle(args.transmission_pkl)
    df_all["global_mean_zscore"] = -df_all["global_mean_zscore"]

    parent = df_filtered[df_filtered["FILTER"] == "parent"]
    grandparent = df_filtered[df_filtered["FILTER"] == "grandparent"]

    all_K_values = sorted(df_all["K"].unique())
    vt_rates = sorted(df_all["PM"].unique())

    pm_df_long = df_all.copy()
    pm_df_long["K"] = pd.Categorical(pm_df_long["K"], categories=all_K_values, ordered=True)
    pm_df_long["PM_label"] = pm_df_long["PM"].map(lambda x: f"{x:.1f}")
    pm_df_long["PM_label"] = pd.Categorical(
        pm_df_long["PM_label"],
        categories=[f"{x:.1f}" for x in vt_rates],
        ordered=True,
    )

    parent_plot = parent.copy()
    parent_plot["Filter"] = "1.0\n(No siblings)"

    grandparent_plot = grandparent.copy()
    grandparent_plot["Filter"] = "1.0\n(No cousins)"

    filter_df_long = pd.concat([parent_plot, grandparent_plot], axis=0)
    filter_df_long["K"] = pd.Categorical(filter_df_long["K"], categories=all_K_values, ordered=True)
    filter_df_long["Filter"] = pd.Categorical(
        filter_df_long["Filter"],
        categories=["1.0\n(No siblings)", "1.0\n(No cousins)"],
        ordered=True,
    )

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

    fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN), dpi=600, layout="constrained")
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 2], figure=fig)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)

    data = all_K_values
    colors = met_brew(name="Homer2", brew_type="continuous", n=6)[3:]

    palette_K = dict(zip(data, colors))

    sns.boxplot(
        data=pm_df_long,
        x="PM_label",
        y="global_mean_zscore",
        hue="K",
        palette=palette_K,
        ax=ax1,
        linewidth=0.75,
        showfliers=True,
        flierprops=dict(
            marker=".",
            markersize=1,
            markerfacecolor="black",
            markeredgecolor="black",
        ),
    )

    ax1.axhline(0, color="black", linestyle="--", linewidth=1)
    ax1.legend_.remove()
    ax1.tick_params(axis="both", which="both", length=0, labelsize=FONT_SIZE)
    ax1.set_xlabel(r"Vertical transmission rate", fontsize=FONT_SIZE)
    ax1.set_ylabel(r"Cophylogeny score", fontsize=FONT_SIZE)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    sns.boxplot(
        data=filter_df_long,
        x="Filter",
        y="global_mean_zscore",
        hue="K",
        palette=palette_K,
        ax=ax2,
        showfliers=True,
        flierprops=dict(
            marker=".",
            markersize=1,
            markerfacecolor="black",
            markeredgecolor="black",
        ),
    )

    ax2.legend(
        loc="upper right",
        bbox_to_anchor=(1.6, 1),
        borderaxespad=0,
        title="Population size",
        fontsize=FONT_SIZE,
        title_fontsize=FONT_SIZE,
        frameon=False,
        edgecolor="none",
        framealpha=0,
    )

    ax2.axhline(0, color="black", linestyle="--", linewidth=1)
    ax2.tick_params(axis="both", which="both", length=0, labelsize=FONT_SIZE)
    ax2.set_xlabel("")
    ax2.set_ylabel("")
    ax2.tick_params(labelleft=False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    fig.savefig(args.output_pdf)
    plt.close(fig)

if __name__ == "__main__":
    main()