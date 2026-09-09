import argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
from scipy.stats import norm

FONT_SIZE = 7
FIG_WIDTH_IN = 150 / 25.4
FIG_HEIGHT_IN = 60 / 25.4

SUPERPOP_MAP = {
    "AFR": ["ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"],
    "AMR": ["CLM", "MXL", "PEL", "PUR"],
    "EAS": ["CDX", "CHB", "CHS", "JPT", "KHV"],
    "EUR": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "SAS": ["BEB", "GIH", "ITU", "PJL", "STU"]
}

def get_superpop(pop):
    for sp, members in SUPERPOP_MAP.items():
        if pop in members:
            return sp
    return None

def get_marker(row):
    if row['p_value'] < 0.001:
        return '***'
    elif row['p_value'] < 0.01:
        return '**'
    elif row['p_value'] < 0.05:
        return '*'
    else:
        return ''

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--whole_genome_zscores",
        required=True,
        help="whole genome zscores file"
    )
    parser.add_argument(
        "--output_figure",
        required=True,
        help="output figure file"
    )
    args = parser.parse_args()

    single_df = pd.read_csv(
        args.whole_genome_zscores,
        sep="\t",
        header=None,
        skiprows=1,
        names=["population", "p_value", "z_score", "lower", "upper"]
    )

    single_df = single_df[~single_df['population'].str.contains("_")].copy()
    single_df['superpop'] = single_df['population'].apply(get_superpop)
    alpha = 0.05 / len(single_df)
    single_df['z_score_upper'] = ((single_df['upper'] - single_df['z_score']) / 1.96) * norm.ppf(1-alpha)
    single_df['z_score_lower'] = ((single_df['z_score'] - single_df['lower']) / 1.96) * norm.ppf(1-alpha)

    base_cmaps = {
        "AFR": cm.Oranges,
        "AMR": cm.Reds,
        "EAS": cm.Greens,
        "EUR": cm.Blues,
        "SAS": cm.Purples
    }

    color_map = {}
    for superpop, pops in SUPERPOP_MAP.items():
        n = len(pops)
        colormap = base_cmaps[superpop]
        colors = [colormap(0.3 + 0.6 * i / (max(n - 1, 1))) for i in range(n)]
        for pop, color in zip(pops, colors):
            color_map[pop] = color

    single_df['color'] = single_df['population'].map(color_map)
    single_df = single_df.sort_values(['superpop', 'z_score'], ascending=[True, False])
    single_df['marker'] = single_df.apply(get_marker, axis=1)

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

    plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN), layout="constrained", dpi=600)
    ax = plt.gca()
    y_positions = np.arange(len(single_df))

    ax.scatter(
        y_positions,
        -single_df['z_score'],
        color=single_df['color'],
        edgecolor='black'
    )

    ax.errorbar(
        y_positions,
        -single_df['z_score'],
        yerr=[single_df['z_score_lower'], single_df['z_score_upper']],
        ecolor='black',
        fmt='none',
        capsize=3,
        zorder=0
    )

    for y, z, marker in zip(y_positions, -single_df['z_score'], single_df['marker']):
        if marker:
            ax.text(
                y - 0.2,
                10.7,
                marker,
                va='bottom',
                ha='left',
                fontsize=FONT_SIZE,
                fontweight='bold',
                color='black'
            )

    ax.set_xticks(y_positions)
    ax.set_xticklabels(single_df['population'], fontsize=FONT_SIZE, rotation=45)
    ax.set_ylabel("Cophylogeny score", fontsize=FONT_SIZE)
    ax.set_xlabel("1000 Genomes Project", fontsize=FONT_SIZE)
    ax.tick_params(axis='x', labelsize=FONT_SIZE)
    ax.tick_params(axis='y', labelsize=FONT_SIZE)
    ax.axhline(0.0, ls='--', color='grey', alpha=0.8, zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.savefig(args.output_figure, format="pdf", bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    main()