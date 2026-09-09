import argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
from scipy.stats import norm

FONT_SIZE = 7
FIG_WIDTH_IN = 140 / 25.4
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--whole_genome_zscores",
        required=True,
        help="whole genome zscores file"
    )
    parser.add_argument(
        "--all_populations_fis",
        required=True,
        help="all populations fis file"
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

    fis_df = pd.read_csv(args.all_populations_fis, sep="\t")
    fis_df.rename(columns={'pop': 'population', 'theta': 'fis'}, inplace=True)
    fis_df = fis_df[~fis_df['population'].str.contains("_")].copy()
    alpha = 0.05 / len(single_df)
    fis_df['fis_upper'] = ((fis_df['upper'] - fis_df['fis']) / 1.96) * norm.ppf(1-alpha)
    fis_df['fis_lower'] = ((fis_df['fis'] - fis_df['lower']) / 1.96) * norm.ppf(1-alpha)
    single_df = pd.merge(single_df, fis_df[['population', 'fis', 'fis_upper', 'fis_lower']], on='population', how='inner')

    single_df['superpop'] = single_df['population'].apply(get_superpop)

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

    ax.errorbar(
        y_positions,
        single_df['fis'],
        yerr=[single_df['fis_lower'], single_df['fis_upper']],
        ecolor='black',
        fmt='none',
        capsize=3,
        zorder=0
    )

    ax.scatter(
        y_positions,
        single_df['fis'],
        color=single_df['color'],
        edgecolor='black'
    )

    ax.set_xticks(y_positions)
    ax.set_xticklabels(single_df['population'], fontsize=FONT_SIZE, rotation=45)
    ax.set_ylabel(r"$F_{IS}$", fontsize=FONT_SIZE)
    ax.set_xlabel("1000 Genomes Project", fontsize=FONT_SIZE)
    ax.tick_params(axis='x', labelsize=FONT_SIZE)
    ax.axhline(0.0, ls='--', color='grey', alpha=0.8, zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.savefig(args.output_figure, format="pdf", bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    main()