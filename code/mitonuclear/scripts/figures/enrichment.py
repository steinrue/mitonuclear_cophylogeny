import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

FONT_SIZE = 7
FIG_WIDTH_IN = 60 / 25.4
FIG_HEIGHT_IN = 100 / 25.4


def pval_to_stars(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return ''


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--enrichment_df_pkl",
        required=True,
        help="enrichment df pkl file",
    )
    parser.add_argument(
        "--output_pdf",
        required=True,
        help="output pdf file",
    )
    args = parser.parse_args()

    enrichment_df = pd.read_pickle(args.enrichment_df_pkl)

    df_1pct = enrichment_df[enrichment_df["thresh"] == 0.99].copy()
    df_1pct["pval"] = df_1pct["pval"] * 26 * 6 # multiple testing correction

    enrich_mat = df_1pct.pivot(index="population", columns="ind", values="enrichment")
    pval_mat = df_1pct.pivot(index="population", columns="ind", values="pval")

    inds_order = [
        "in_exon",
        "in_mitocarta_exon",
        "in_oxphos_exon",
        "in_gene",
        "in_mitocarta_gene",
        "in_oxphos_gene",
    ]
    enrich_mat = enrich_mat[inds_order]
    pval_mat = pval_mat[inds_order]

    annot_mat = pval_mat.map(pval_to_stars)

    log2_enrich_mat = np.log2(enrich_mat)

    single_pops_order = [
        'YRI', 'ESN', 'MSL', 'GWD', 'ACB', 'LWK', 'ASW',
        'PEL', 'PUR', 'CLM', 'MXL',
        'JPT', 'CHS', 'KHV', 'CHB', 'CDX',
        'TSI', 'GBR', 'FIN', 'IBS', 'CEU',
        'BEB', 'GIH', 'PJL', 'STU', 'ITU'
    ]

    single_pops_order = [p for p in single_pops_order if p in log2_enrich_mat.index]
    single_pops_order = single_pops_order[::-1]

    log2_enrich_single = log2_enrich_mat.reindex(single_pops_order)
    annot_single = annot_mat.reindex(single_pops_order)

    vmin = -2.0
    vmax = 2.0
    cmap = "bwr"

    xtick_labels = [
        "Exons",
        "MitoCarta exons",
        "OXPHOS exons",
        "Genes",
        "MitoCarta genes",
        "OXPHOS genes",
    ]

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

    fig1, ax1 = plt.subplots(
        1,
        1,
        figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
        layout="constrained",
        dpi=600,
    )

    h1 = sns.heatmap(
        log2_enrich_single,
        ax=ax1,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        annot=annot_single,
        fmt="",
        cbar=False,
        linewidths=0.5,
        linecolor="black"
    )

    ax1.set_title("Enrichment in top 1% cophylogeny", fontsize=FONT_SIZE)
    ax1.set_xlabel("", fontsize=FONT_SIZE)
    ax1.set_ylabel("1000 Genomes Project", fontsize=FONT_SIZE)

    ax1.set_xticklabels(xtick_labels, rotation=45, ha="right", fontsize=FONT_SIZE)
    ax1.tick_params(axis='y', rotation=0, labelsize=FONT_SIZE)

    cbar1 = fig1.colorbar(h1.collections[0], ax=ax1, shrink=0.4)
    cbar1.set_label("log2(enrichment ratio)", fontsize=FONT_SIZE)
    cbar1.ax.tick_params(labelsize=FONT_SIZE)

    plt.savefig(args.output_pdf, format="pdf", bbox_inches="tight", dpi=600)
    plt.close()


if __name__ == "__main__":
    main()