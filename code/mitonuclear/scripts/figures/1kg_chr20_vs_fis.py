import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress
from scipy.stats import norm

FONT_SIZE = 7
FIG_WIDTH_IN = 100 / 25.4
FIG_HEIGHT_IN = 80 / 25.4

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--all_populations_fis",
        required=True,
        help="all populations fis file"
    )
    parser.add_argument(
        "--chr20_zscores",
        required=True,
        help="chr20 zscores file"
    )
    parser.add_argument(
        "--output_figure",
        required=True,
        help="output figure file"
    )
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

    plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN), layout="constrained", dpi=600)
    ax = plt.gca()
    
    results_df = pd.read_csv(args.chr20_zscores, sep="\t", header=0)
    chr20_df = results_df[results_df['population'].str.contains("_")].copy()
    fis_df = pd.read_csv(args.all_populations_fis, sep="\t")
    alpha = 0.05 / len(chr20_df)
    fis_df['fis_upper'] = ((fis_df['upper'] - fis_df['fis']) / 1.96) * norm.ppf(1-alpha)
    fis_df['fis_lower'] = ((fis_df['fis'] - fis_df['lower']) / 1.96) * norm.ppf(1-alpha)
    chr20_df['z_score_upper'] = ((chr20_df['upper'] - chr20_df['z_score']) / 1.96) * norm.ppf(1-alpha)
    chr20_df['z_score_lower'] = ((chr20_df['z_score'] - chr20_df['lower']) / 1.96) * norm.ppf(1-alpha)

    dual_df = pd.merge(
        chr20_df,
        fis_df[['population', 'fis', 'fis_lower', 'fis_upper']],
        on='population',
        how='inner'
    )

    paired_all_chrs = ["YRI_TSI", "IBS_TSI", "GIH_STU", "CHS_GIH", "JPT_IBS"]
    sampled = dual_df[dual_df['population'].isin(paired_all_chrs)].copy()
    
    scatter = ax.scatter(
        dual_df['fis'],
        -dual_df['z_score'],
        color='slategrey',
        alpha=1,
        s=10,
        edgecolors='none',
        label="Combined populations"
    )

    ax.errorbar(
        dual_df['fis'],
        -dual_df['z_score'],
        xerr=[dual_df['fis_lower'], dual_df['fis_upper']],
        yerr=[dual_df['z_score_lower'], dual_df['z_score_upper']],
        ecolor='slategrey',
        fmt='none',
        linewidth=0.6,
        capsize=0.3,
        zorder=0
    )

    for xi, yi, label in zip(sampled['fis'], -sampled['z_score'], sampled['population']):
            display_label = label.replace('_', ' and ')
            ax.annotate(
                display_label,
                xy=(xi, yi),                 
                xytext=(xi + 0.005, yi - 5),  
                va='bottom',
                ha='left',
                fontsize=FONT_SIZE,
                color='black',
                arrowprops=dict(
                    arrowstyle="-",
                    color="black",
                    linewidth=0.5
                ))

    slope, intercept, r, p, _ = linregress(dual_df['fis'], -dual_df['z_score'])
    x_line = np.linspace(dual_df['fis'].min(), 0.1, 100)
    y_line = intercept + slope * x_line
    reg_line, = ax.plot(
        x_line,
        y_line,
        linestyle='--',
        color='red',
        lw=1.2,
        label=f'Regression: $R^2$={r**2:.2f}'
    )

    ax.legend(
        handles=[scatter, reg_line],
        loc='upper left',
        frameon=True,
        fontsize=FONT_SIZE
    )

    ax.set_xlabel(r"$F_{IS}$", fontsize=FONT_SIZE)
    ax.set_ylabel("Cophylogeny score", fontsize=FONT_SIZE)
    ax.tick_params(axis='x', labelsize=FONT_SIZE)
    ax.tick_params(axis='y', labelsize=FONT_SIZE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.savefig(args.output_figure, format="pdf", bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    main()