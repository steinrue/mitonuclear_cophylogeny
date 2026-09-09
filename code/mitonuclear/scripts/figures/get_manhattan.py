import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

FONT_SIZE = 7
FIG_WIDTH_IN = 140 / 25.4
FIG_HEIGHT_IN = 35 / 25.4

SUPERPOP_MAP = {
    "AFR": ["ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"],
    "AMR": ["CLM", "MXL", "PEL", "PUR"],
    "EAS": ["CDX", "CHB", "CHS", "JPT", "KHV"],
    "EUR": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "SAS": ["BEB", "GIH", "ITU", "PJL", "STU"]
}

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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--single_locus_pvals",
        required=True,
        help="single locus pvals file"
    )
    parser.add_argument(
        "--output_figure",
        required=True,
        help="output figure file"
    )
    args = parser.parse_args()

    population = os.path.basename(args.single_locus_pvals).split(".")[0]

    df_pvals = pd.read_csv(args.single_locus_pvals, sep="\t")
    df_pvals['mid_pos'] = (df_pvals['start_pos'] + df_pvals['end_pos']) / 2
    ordered_chrs = [f"chr{i}" for i in range(1, 23)]
    df_pvals = df_pvals[df_pvals['chr'].isin(ordered_chrs)].copy()

    color = color_map[population]
    light_color = tuple(0.6 * c + 0.4 for c in color[:3])

    chr_to_order = {chrom: i for i, chrom in enumerate(ordered_chrs)}
    df_pvals['chr_order'] = df_pvals['chr'].map(chr_to_order)

    chr_offsets = {}
    current_offset = 0
    for chrom in ordered_chrs:
        chr_data = df_pvals[df_pvals['chr'] == chrom]
        chr_offsets[chrom] = current_offset
        chr_len = chr_data['mid_pos'].max()
        current_offset += chr_len

    df_pvals['offset'] = df_pvals['chr'].map(chr_offsets)
    df_pvals['cum_pos'] = df_pvals['mid_pos'] + df_pvals['offset']

    grey_color_map = {}
    for chrom in ordered_chrs:
        idx = chr_to_order[chrom]
        grey_color_map[chrom] = light_color if idx % 2 == 0 else color
    df_pvals['color'] = df_pvals['chr'].map(grey_color_map)

    genome_wide_thr = -np.log10(0.05 / len(df_pvals))

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
    ax.scatter(
        df_pvals['cum_pos'],
        df_pvals['log10_pval'],
        c=df_pvals['color'],
        s=3,
        rasterized=True,
        linewidths=0
    )

    ax.axhline(
        genome_wide_thr,
        color='grey',
        linestyle='--',
        linewidth=0.7,
        label="BH single population"
    )

    ax.axhline(
        -np.log10(4.8e-9),
        color='black',
        linestyle='--',
        linewidth=0.7,
        label="BH multiple populations"
    )

    chr_ticks = []
    chr_labels = []
    for chrom in ordered_chrs:
        chr_data = df_pvals[df_pvals['chr'] == chrom]
        center = (chr_data['cum_pos'].min() + chr_data['cum_pos'].max()) / 2
        chr_ticks.append(center)
        chr_labels.append(chrom)

    ax.set_xticks(chr_ticks)
    ax.set_xticklabels(chr_labels, rotation=90, fontsize=FONT_SIZE)

    ax.set_ylabel(r"$-\log_{10}(p)$", fontsize=FONT_SIZE)
    ax.set_xlabel("Chromosome", fontsize=FONT_SIZE)
    ax.tick_params(axis='y', labelsize=FONT_SIZE)
    ax.tick_params(axis='x', labelsize=FONT_SIZE)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    #ax.legend(fontsize=FONT_SIZE, loc='upper right', frameon=False, bbox_to_anchor=(1, 1.15))
    plt.ylim(0, 9.5)
    plt.title(f"Mitonuclear cophylogeny in {population}", fontsize=FONT_SIZE)

    plt.savefig(args.output_figure, format="pdf", bbox_inches="tight", dpi = 600)
    plt.close()

if __name__ == "__main__":
    main()