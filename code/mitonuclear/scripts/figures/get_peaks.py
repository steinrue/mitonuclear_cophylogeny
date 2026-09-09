import argparse
import gzip
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np

FONT_SIZE = 7
FIG_WIDTH_IN = 100 / 25.4
FIG_HEIGHT_IN = 45 / 25.4

CHR_LENGTHS_GRCH38 = {
    "chr1":   248956422,
    "chr2":   242193529,
    "chr3":   198295559,
    "chr4":   190214555,
    "chr5":   181538259,
    "chr6":   170805979,
    "chr7":   159345973,
    "chr8":   145138636,
    "chr9":   138394717,
    "chr10":  133797422,
    "chr11":  135086622,
    "chr12":  133275309,
    "chr13":  114364328,
    "chr14":  107043718,
    "chr15":  101991189,
    "chr16":   90338345,
    "chr17":   83257441,
    "chr18":   80373285,
    "chr19":   58617616,
    "chr20":   64444167,
    "chr21":   46709983,
    "chr22":   50818468
}

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

color_map_1kg = {}
for superpop, pops in SUPERPOP_MAP.items():
    n = len(pops)
    colormap = base_cmaps[superpop]
    colors = [colormap(0.3 + 0.6 * i / (max(n - 1, 1))) for i in range(n)]
    for pop, color in zip(pops, colors):
        color_map_1kg[pop] = color


def parse_gff3_protein_coding_df(gff3_path):
    rows = []
    with gzip.open(gff3_path, 'rt') as gff:
        for line in gff:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom = fields[0]
            feature_type = fields[2]
            if chrom not in CHR_LENGTHS_GRCH38:
                continue
            if feature_type not in ('gene', 'exon'):
                continue
            start_1based = fields[3]
            end_1based = fields[4]
            attr = fields[8]
            gene_id = None
            gene_name = None
            gene_type = None
            for kv in attr.split(';'):
                key, value = kv.split('=', 1)
                if key == 'gene_id':
                    gene_id = value.split('.')[0]
                elif key == 'gene_name':
                    gene_name = value
                elif key == 'gene_type':
                    gene_type = value
            start_0based = int(start_1based) - 1
            end_0based = int(end_1based)
            rows.append({
                'chr': chrom,
                'feature': feature_type,
                'start': start_0based,
                'end': end_0based,
                'gene_id': gene_id,
                'gene_name': gene_name,
                'gene_type': gene_type
            })
    df = pd.DataFrame(
        rows,
        columns=['chr', 'feature', 'start', 'end', 'gene_id', 'gene_name', 'gene_type']
    )
    df = df[df['gene_type'] == 'protein_coding']
    return df


def get_superpop(pop):
    for sp, members in SUPERPOP_MAP.items():
        if pop in members:
            return sp
    return None


def _get_peak_region(df_pvals, flank_bp):
    peak_row = df_pvals.loc[df_pvals['log10_pval'].idxmax()]
    chr_peak = peak_row['chr']
    peak_mid = 0.5 * (peak_row['start_pos'] + peak_row['end_pos'])
    region_start = peak_mid - flank_bp
    region_end = peak_mid + flank_bp
    return chr_peak, peak_mid, region_start, region_end


def plot_peak_region_multi_pop(
    df_pvals_dict,
    gene_df,
    populations,
    mito_genes,
    oxphos_genes,
    flank_bp=100_000,
    fig_width_in=FIG_WIDTH_IN,
    fig_height_in=FIG_HEIGHT_IN,
    dpi=600
):
    n_pops = len(populations)

    fig, axes = plt.subplots(
        1,
        n_pops,
        figsize=(fig_width_in, fig_height_in),
        dpi=dpi,
        sharey=True,
        layout="constrained"
    )

    if n_pops == 1:
        axes = [axes]

    color_gw = '0.5'
    color_sw = 'black'

    legend_handles = []
    legend_labels = []

    for ax, pop in zip(axes, populations):
        df_pvals = df_pvals_dict[pop]

        if isinstance(flank_bp, dict):
            flank_this = flank_bp.get(pop, 100_000)
        else:
            flank_this = flank_bp

        chr_peak, peak_mid, region_start, region_end = _get_peak_region(df_pvals, flank_this)

        df_region = df_pvals[
            (df_pvals['chr'] == chr_peak) &
            (df_pvals['start_pos'] <= region_end) &
            (df_pvals['end_pos'] >= region_start)
        ].copy()

        if df_region.empty:
            raise ValueError(f"No pvals found in the specified region around the peak for {pop}.")

        df_region['mid_pos'] = 0.5 * (df_region['start_pos'] + df_region['end_pos'])

        df_chr = df_pvals[df_pvals['chr'] == chr_peak]
        chr_mean_log10_pval = df_chr['log10_pval'].mean()

        N_tests = len(df_pvals)
        genome_wide_threshold = -np.log10(0.05 / N_tests)
        study_wide_threshold = -np.log10(4.8e-9)

        color_pop = color_map_1kg.get(pop, "black")
        main_line = ax.plot(
            df_region['mid_pos'] / 1e6,
            df_region['log10_pval'],
            color=color_pop,
            linewidth=1.0
        )[0]

        hl_gw = ax.axhline(
            genome_wide_threshold,
            color=color_gw,
            linestyle='--',
            linewidth=0.3
        )
        hl_sw = ax.axhline(
            study_wide_threshold,
            color=color_sw,
            linestyle='--',
            linewidth=0.3
        )

        genes_region = gene_df[
            (gene_df['chr'] == chr_peak) &
            (gene_df['end'] >= region_start) &
            (gene_df['start'] <= region_end)
        ].copy()

        ax.autoscale()
        ylim_low, ylim_high = ax.get_ylim()

        for _, gene in genes_region.iterrows():
            g_start = max(gene['start'], region_start)
            g_end = min(gene['end'], region_end)

            base_patch = ax.axvspan(
                g_start / 1e6,
                g_end / 1e6,
                color='0.8',
                alpha=0.3,
                linewidth=0
            )

            if gene['gene_id'] in mito_genes:
                ax.axvspan(
                    g_start / 1e6,
                    g_end / 1e6,
                    color='orange',
                    alpha=0.3,
                    linewidth=0
                )

            if gene['gene_id'] in oxphos_genes:
                ax.axvspan(
                    g_start / 1e6,
                    g_end / 1e6,
                    color='red',
                    alpha=0.3,
                    linewidth=0
                )

                x_pos = 0.5 * (g_start + g_end) / 1e6
                y_pos = ylim_low
                ax.text(
                    x_pos,
                    y_pos,
                    gene['gene_name'] if pd.notnull(gene['gene_name']) else gene['gene_id'],
                    ha='center',
                    va='center',
                    fontsize=FONT_SIZE,
                    rotation=0,
                    clip_on=True
                )

        ax.set_xlim(region_start / 1e6, region_end / 1e6)
        ax.set_ylim(-0.7, 9)

        ax.set_xlabel(f"{chr_peak} position (Mb)", fontsize=FONT_SIZE)
        ax.set_title(pop, fontsize=FONT_SIZE)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.tick_params(axis='both', which='both', length=2, width=0.3, direction='out', labelsize=FONT_SIZE)

        if not legend_handles:
            legend_handles.append(main_line)
            legend_labels.append("Mitonuclear cophylogeny")

            legend_handles.extend([hl_gw, hl_sw])
            legend_labels.extend([
                "Genome-wide significance\n(pop.-specific)",
                "Genome-wide significance\n(across pops.)"
            ])

            all_genes_patch = mpatches.Patch(color='0.8', alpha=0.5)
            mito_patch = mpatches.Patch(color='orange', alpha=0.5)
            oxphos_patch = mpatches.Patch(color='red', alpha=0.5)

            legend_handles.extend([all_genes_patch, mito_patch, oxphos_patch])
            legend_labels.extend([
                "All genes",
                "MitoCarta genes",
                "OXPHOS genes"
            ])

    axes[0].set_ylabel(r"$\log_{10}(\mathrm{p\ value})$", fontsize=FONT_SIZE)

    fig.legend(
        legend_handles,
        legend_labels,
        fontsize=FONT_SIZE,
        frameon=False,
        loc='upper right',
        bbox_to_anchor=(1.4, 0.8)
    )

    return fig, axes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--whole_genome_zscores",
        required=True,
        help="whole genome zscores file"
    )
    parser.add_argument(
        "--gencode_gff3_gz",
        required=True,
        help="gencode gff3 gz file"
    )
    parser.add_argument(
        "--mitocarta_xls",
        required=True,
        help="mitocarta xls file"
    )
    parser.add_argument(
        "--jpt_pvals",
        required=True,
    )
    parser.add_argument(
        "--ibs_pvals",
        required=True,
    )
    parser.add_argument(
        "--output_pdf",
        required=True,
        help="output pdf file"
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

    single_df = pd.read_csv(args.whole_genome_zscores, sep="\t")
    all_populations = single_df["population"].to_list()
    single_df = single_df[(~single_df['population'].str.contains("_"))].copy()

    df = parse_gff3_protein_coding_df(args.gencode_gff3_gz)
    gene_df = df[df['feature'] == 'gene']

    mito_df = pd.read_excel(args.mitocarta_xls, sheet_name=1)
    mito_genes = set()
    for ids in mito_df['EnsemblGeneID_mapping_version_20200130'].dropna().astype(str):
        for gene_id in ids.split('|'):
            mito_genes.add(gene_id.split('.')[0])

    oxphos_genes = set()
    oxphos_rows = mito_df.loc[
        mito_df["MitoCarta3.0_MitoPathways"].str.contains("OXPHOS", na=False, regex=False),
        "EnsemblGeneID_mapping_version_20200130"
    ].dropna().astype(str)
    for ids in oxphos_rows:
        for gene_id in ids.split('|'):
            oxphos_genes.add(gene_id.split('.')[0])

    df_pvals_dict = {}
    df_pvals_dict["JPT"] =  pd.read_csv(args.jpt_pvals, sep="\t")
    df_pvals_dict["IBS"] =  pd.read_csv(args.ibs_pvals, sep="\t")

    populations = ["JPT", "IBS"]
    flanks = {"JPT": 270_000, "IBS": 550_000}

    fig, axes = plot_peak_region_multi_pop(
        df_pvals_dict=df_pvals_dict,
        gene_df=gene_df,
        populations=populations,
        mito_genes=mito_genes,
        oxphos_genes=oxphos_genes,
        flank_bp=flanks
    )
    plt.savefig(args.output_pdf, format="pdf", bbox_inches="tight", dpi=600)
    plt.close()


if __name__ == "__main__":
    main()