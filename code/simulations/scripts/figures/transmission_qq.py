import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import beta
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from met_brewer import met_brew

FONT_SIZE = 7
FIG_WIDTH_IN = 100 / 25.4
FIG_HEIGHT_IN = 50 / 25.4

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--transmission-pkl", required=True)
    parser.add_argument("--transmission-w-filter-pkl", required=True)
    parser.add_argument("--output-pdf", required=True)
    args = parser.parse_args()

    np.random.seed(42)

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
        1, 2,
        figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
        sharey=True,
        layout="constrained",
        dpi=600,
    )

    sample_step = 10000
    colors = met_brew(name="Homer2", brew_type="continuous", n=6)[3:]

    df1 = pd.read_pickle(args.transmission_pkl)
    vertical = df1[df1["PM"] == 1.0]
    grouped = vertical.groupby("K")
    data, labels = [], []
    for m_val, group in grouped:
        sampled_vectors = []
        for row in group.itertuples():
            pval_list = [float(x) for x in row.p_vals.split(";")]
            start_list = [int(float(x)) for x in row.start_vals.split(";")]
            idxs = np.array([(np.searchsorted(start_list, pos, side="right") - 1) for pos in range(0, int(1e7), sample_step)])
            sampled_vectors.append(np.sort([pval_list[idx] for idx in idxs]))
        data.append(np.median(np.vstack(sampled_vectors), axis=0))
        labels.append(f"{m_val}")
    
    n = len(data[0])
    exp = -np.log10((n - np.arange(1, n + 1) + 1) / n)
    i, m = np.arange(1, n + 2), n + 2
    lower, upper = -np.log10(beta.ppf(0.025, i, m - i)), -np.log10(beta.ppf(0.975, i, m - i))
    axs[0].fill_between(-np.log10(i / m), lower, upper, color="k", alpha=0.1)
    for vals, label, color in zip(data, labels, colors):
        obs = np.sort(vals)
        axs[0].scatter(exp[:len(obs)], obs, label=label, alpha=0.8, s=8, color=color, edgecolor="none", rasterized=True)
    axs[0].plot([0, 5], [0, 5], linestyle="--", color="grey", linewidth=1)
    axs[0].set_xlim(0, 3.5)
    axs[0].set_ylim(0, 5)
    axs[0].spines["top"].set_alpha(0)
    axs[0].spines["right"].set_alpha(0)
    axs[0].set_ylabel(r"$-\log_{10}(p)$", fontsize=FONT_SIZE)
    axs[0].set_title(r"No filter", fontsize=FONT_SIZE)
    #axs[0].legend(loc="upper left", borderaxespad=0, title="Population size", fontsize=FONT_SIZE, title_fontsize=FONT_SIZE, frameon=False, edgecolor="none", framealpha=0)
    axs[0].tick_params(axis="both", labelsize=FONT_SIZE, length=2)
    axs[0].grid(False)
    axins0 = inset_axes(axs[0], width="30%", height="30%", loc="lower right", bbox_to_anchor=(0.0, 0.05, 0.92, 0.92), bbox_transform=axs[0].transAxes, borderpad=0)
    for vals, color in zip(data, colors):
        axins0.scatter(np.sort(10 ** (-exp))[:len(np.sort(10 ** (-np.array(vals))))], np.sort(10 ** (-np.array(vals))), s=3, alpha=0.8, color=color, edgecolor="none", rasterized=True)
    axins0.plot([0, 1], [0, 1], linestyle="--", color="grey", linewidth=1)
    axins0.set_xlim(0, 1)
    axins0.set_ylim(0, 1)
    axins0.tick_params(axis="y", which="both", left=False)
    axins0.yaxis.tick_right()
    axins0.yaxis.set_label_position("right")
    axins0.set_xticks([])
    axins0.tick_params(axis="both", labelsize=FONT_SIZE, length=2)
    axins0.grid(False)

    df2 = pd.read_pickle(args.transmission_w_filter_pkl)
    vertical = df2[df2["FILTER"] == "parent"]
    grouped = vertical.groupby("K")
    data, labels = [], []
    for m_val, group in grouped:
        sampled_vectors = []
        for row in group.itertuples():
            pval_list = [float(x) for x in row.p_vals.split(";")]
            start_list = [int(float(x)) for x in row.start_vals.split(";")]
            idxs = np.array([(np.searchsorted(start_list, pos, side="right") - 1) for pos in range(0, int(1e7), sample_step)])
            sampled_vectors.append(np.sort([pval_list[idx] for idx in idxs]))
        data.append(np.median(np.vstack(sampled_vectors), axis=0))
        labels.append(f"{m_val}")
    
    axs[1].fill_between(-np.log10(i / m), lower, upper, color="k", alpha=0.1)
    for vals, label, color in zip(data, labels, colors):
        obs = np.sort(vals)
        axs[1].scatter(exp[:len(obs)], obs, label=label, alpha=0.8, s=8, color=color, edgecolor="none", rasterized=True)
    axs[1].plot([0, 5], [0, 5], linestyle="--", color="grey", linewidth=1)
    axs[1].set_xlim(0, 3.5)
    axs[1].spines["top"].set_alpha(0)
    axs[1].spines["right"].set_alpha(0)
    axs[1].set_title(r"No siblings", fontsize=FONT_SIZE)
    #axs[1].legend(loc="upper left", borderaxespad=0, title="Population size", fontsize=FONT_SIZE, title_fontsize=FONT_SIZE, frameon=False, edgecolor="none", framealpha=0)
    axs[1].tick_params(axis="both", labelsize=FONT_SIZE, length=2)
    axs[1].grid(False)
    axins1 = inset_axes(axs[1], width="30%", height="30%", loc="lower right", bbox_to_anchor=(0.0, 0.05, 0.92, 0.92), bbox_transform=axs[1].transAxes, borderpad=0)
    for vals, color in zip(data, colors):
        axins1.scatter(np.sort(10 ** (-exp))[:len(np.sort(10 ** (-np.array(vals))))], np.sort(10 ** (-np.array(vals))), s=3, alpha=0.8, color=color, edgecolor="none", rasterized=True)
    axins1.plot([0, 1], [0, 1], linestyle="--", color="grey", linewidth=1)
    axins1.set_xlim(0, 1)
    axins1.set_ylim(0, 1)
    axins1.tick_params(axis="y", which="both", left=False)
    axins1.yaxis.tick_right()
    axins1.yaxis.set_label_position("right")
    axins1.set_xticks([])
    axins1.tick_params(axis="both", labelsize=FONT_SIZE, length=2)
    axins1.grid(False)

    fig.supxlabel(r"Expected $-\log_{10}(p)$", fontsize=FONT_SIZE)

    fig.savefig(args.output_pdf, format="pdf", dpi=600)
    plt.close(fig)

if __name__ == "__main__":
    main()