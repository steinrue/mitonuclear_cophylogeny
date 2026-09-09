import argparse
import pandas as pd
import numpy as np
from cyvcf2 import VCF
from rpy2.robjects import r
from rpy2.robjects.vectors import FloatVector
from rpy2.robjects.packages import importr
import tskit
from ete3 import Tree
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import matplotlib.colors as mcolors
import statsmodels.api as sm

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gencode_path", type=str, required=True, help="gencode file")
    parser.add_argument("--results_df", type=str, required=True, help="results file")
    parser.add_argument("--vcf", type=str, required=True, help="vcf file")
    parser.add_argument("--hsd", type=str, required=True, help="hsd file")
    parser.add_argument("--trees", type=str, required=True, help="trees file")
    parser.add_argument("--samples_list", type=str, required=True, help="samples list file")
    parser.add_argument("--mt_tree", type=str, required=True, help="mt tree file")
    parser.add_argument("--output_pdf", type=str, required=True, help="output pdf file")
    args = parser.parse_args()

    gencode_path = args.gencode_path
    results_df = pd.read_csv(args.results_df, sep="\t")

    chrom = 'chr19'
    region_start = int(results_df[results_df["log10_pval"] == results_df["log10_pval"].max()]["start_pos"])
    region_end = int(results_df[results_df["log10_pval"] == results_df["log10_pval"].max()]["end_pos"])

    vcf = VCF(args.vcf)
    variant_data = []
    variant_ids = []

    for variant in vcf(f'{chrom}:{region_start}-{region_end}'):
        gt_array = np.array(variant.gt_types)
        print(gt_array)
        genotypes = gt_array[gt_array != 2]
        unique_genotypes = np.unique(genotypes)
        if len(unique_genotypes) > 1 and np.any(np.isin(unique_genotypes, [0, 1, 3])):
            variant_ids.append(f'{chrom}_{variant.POS}_{variant.REF}_{variant.ALT[0]}')
            variant_data.append(gt_array)
    gt_df = pd.DataFrame(variant_data, columns=vcf.samples, index=variant_ids)

    hsd_lines = []
    with open(args.hsd, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) >= 3:
                hsd_lines.append({'SampleID': fields[0], 'Haplogroup': fields[2]})

    hsd_df = pd.DataFrame(hsd_lines)
    hsd_df = hsd_df.drop_duplicates('SampleID').set_index('SampleID')
    hsd_df['Haplogroup'] = hsd_df['Haplogroup'].str.extract(r'([A-Z][0-9]*)')

    shared_samples = hsd_df.index.intersection(gt_df.columns)
    gt_df_aligned = gt_df[shared_samples].transpose()
    gt_df_aligned.index.name = 'SampleID'

    haplo_info = hsd_df.loc[shared_samples, ['Haplogroup']]
    haplo_info['SampleID'] = haplo_info.index
    final_df = haplo_info.join(gt_df_aligned, how='inner')
    cols = ['SampleID', 'Haplogroup'] + list(gt_df_aligned.columns)
    final_df = final_df[cols]

    ape = importr('ape')
    TreeDist = importr('TreeDist')
    ts = tskit.load(args.trees)
    with open(args.samples_list, 'r') as f:
        idsToExtract = [line.strip() for line in f]
    tree = ts.at(.5 * (region_start + region_end))
    num_samples = ts.num_samples
    num_individuals = num_samples // 2
    avg_tmrca_matrix = np.zeros((num_individuals, num_individuals))
    for i in range(num_individuals):
        for j in range(i + 1, num_individuals):
            s1, s2 = 2 * i, 2 * i + 1
            t1, t2 = 2 * j, 2 * j + 1
            avg_tmrca_matrix[i, j] = (tree.tmrca(s1, t1) + tree.tmrca(s1, t2) + tree.tmrca(s2, t1) + tree.tmrca(s2, t2)) / 2
            avg_tmrca_matrix[j, i] = avg_tmrca_matrix[i, j]

    r_mat = r['matrix'](FloatVector(avg_tmrca_matrix.flatten()), nrow=num_individuals, ncol=num_individuals)
    r_df = r['as.data.frame'](r_mat)
    r_df.colnames = idsToExtract
    r_df.rownames = idsToExtract
    r.assign("D", r_df)
    r('host_tree <- phangorn::upgma(D)')
    r(f'mt_tree <- read.tree("{args.mt_tree}")')
    r('match_data <- MutualClusteringInfo(host_tree, mt_tree, reportMatching = TRUE)')
    r('split_text <- attr(match_data, "matchedSplits")')
    r('scores <- attr(match_data, "matchedScores")')
    r('top_index <- order(scores, decreasing = TRUE)[1]')
    top_split = str(r('split_text[top_index]'))
    clade_tips = top_split.split("=>")[1].split("|", 1)[0].split()

    tree_path = args.mt_tree
    t = Tree(tree_path, format=1)
    info = final_df.set_index("SampleID")

    if "Haplogroup" not in info.columns:
        raise ValueError("Haplogroup column not found")

    d4_tips = info.index[info["Haplogroup"] == "D4"].tolist()

    clade_tips_present = [leaf.name for leaf in t.iter_leaves() if leaf.name in clade_tips]
    if len(clade_tips_present) == 0:
        raise ValueError("None of the specified clade tips are present in the tree.")

    mrca = t.get_common_ancestor(clade_tips_present)
    mrca_leaf_names = {leaf.name for leaf in mrca.iter_leaves()}
    missing_in_mrca = [tip for tip in clade_tips_present if tip not in mrca_leaf_names]
    if missing_in_mrca:
        raise ValueError("The specified clade does not form a single clade in the tree.")

    all_leaves = list(t.iter_leaves())
    tip_order = [leaf.name for leaf in all_leaves]
    y_positions = {tip: i for i, tip in enumerate(tip_order)}

    def set_x_raw(node, current_dist=0.0):
        node.add_feature("x_raw", current_dist)
        for child in node.children:
            length = child.dist if child.dist is not None else 0.0
            set_x_raw(child, current_dist + length)

    set_x_raw(t)
    max_dist = max(leaf.x_raw for leaf in t.iter_leaves())

    def set_x_final(node):
        if node.is_leaf():
            node.add_feature("x", max_dist)
        else:
            node.add_feature("x", node.x_raw)
            for child in node.children:
                set_x_final(child)

    set_x_final(t)

    def set_y(node):
        if node.is_leaf():
            node.add_feature("y", y_positions[node.name])
        else:
            for child in node.children:
                set_y(child)
            ys = [child.y for child in node.children]
            node.add_feature("y", sum(ys) / len(ys))

    set_y(t)

    clade_leaf_ys = [leaf.y for leaf in mrca.iter_leaves()]
    clade_y_min = min(clade_leaf_ys)
    clade_y_max = max(clade_leaf_ys)
    box_y_min = clade_y_min - 0.5
    box_y_max = clade_y_max + 0.5
    box_x_min = mrca.x
    box_x_max = max_dist

    missing_tips = [tip for tip in tip_order if tip not in info.index]
    if missing_tips:
        extra = pd.DataFrame(index=missing_tips, columns=info.columns)
        info_full = pd.concat([info, extra], axis=0)
    else:
        info_full = info.copy()

    info_full = info_full.loc[tip_order]
    drop_cols = []
    info_full = info_full.drop(columns=[c for c in drop_cols if c in info_full.columns], errors="ignore")

    haplogroups = sorted(info_full["Haplogroup"].dropna().unique())
    n_haplogroups = len(haplogroups)

    palette = plt.get_cmap("tab20")
    hap_colors = {hap: palette(i % 20) for i, hap in enumerate(haplogroups)}

    variant_cols = [c for c in info_full.columns if c != "Haplogroup"]
    variants_numeric = info_full[variant_cols].apply(pd.to_numeric, errors="coerce")

    mat = pd.DataFrame(index=info_full.index)
    mat["Haplogroup"] = info_full["Haplogroup"].values
    for col in variant_cols:
        mat[col] = variants_numeric[col].values

    n_rows, n_cols = mat.shape

    color_matrix = np.empty((n_rows, n_cols), dtype=object)
    for i, tip in enumerate(mat.index):
        for j, var in enumerate(mat.columns):
            if var == "Haplogroup":
                hap = mat.iloc[i, j]
                color_matrix[i, j] = hap_colors.get(hap, "white")
            else:
                val = mat.iloc[i, j]
                if val == 0:
                    color_matrix[i, j] = "white"
                elif val == 1:
                    color_matrix[i, j] = "lightgrey"
                elif val == 3:
                    color_matrix[i, j] = "dimgray"
                else:
                    color_matrix[i, j] = "white"

    rgb_array = np.zeros((n_rows, n_cols, 3))
    for i in range(n_rows):
        for j in range(n_cols):
            rgb_array[i, j, :] = mcolors.to_rgb(color_matrix[i, j])

    info_reg = info_full.copy()
    info_reg['is_D4'] = (info_reg['Haplogroup'] == 'D4').astype(int)
    reg_variant_cols = [c for c in info_reg.columns if c.startswith('chr')]
    info_reg[reg_variant_cols] = info_reg[reg_variant_cols].replace(3, 2)
    y = info_reg['is_D4'].values
    sig_variants = []
    for var in reg_variant_cols:
        x = info_reg[var].values.astype(float)
        X_design = sm.add_constant(x)
        model = sm.OLS(y, X_design).fit()
        slope_p = model.pvalues[1]
        if slope_p < 0.05 / 8:
            sig_variants.append(var)

    FONT_SIZE = 7
    FIG_WIDTH_IN = 80 / 25.4
    FIG_HEIGHT_IN = 180 / 25.4

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

    fig, (ax_tree, ax_heat) = plt.subplots(
        ncols=2,
        figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
        layout="constrained",
        dpi=600,
        sharey=True
    )

    rect = Rectangle(
        (box_x_min, box_y_min),
        box_x_max - box_x_min,
        box_y_max - box_y_min,
        facecolor="lightgrey",
        edgecolor="none",
        zorder=0
    )
    ax_tree.add_patch(rect)

    for node in t.traverse():
        if not node.is_root():
            parent = node.up
            ax_tree.plot([parent.x, parent.x], [parent.y, node.y],
                         color="black", linewidth=0.5, zorder=1)
            ax_tree.plot([parent.x, node.x], [node.y, node.y],
                         color="black", linewidth=0.5, zorder=1)

    for leaf in t.iter_leaves():
        color = "black"
        ax_tree.text(leaf.x, leaf.y, f" {leaf.name}",
                     va="center", ha="left",
                     color=color, fontsize=5, zorder=2)

    ax_heat.imshow(rgb_array, aspect="auto", origin="upper")
    
    ax_heat.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax_heat.grid(which="minor", color="black", linestyle="-", linewidth=0.5)
    ax_heat.tick_params(which="minor", bottom=False, left=False)

    ax_heat.set_yticks(np.arange(n_rows))
    ax_heat.set_yticklabels(mat.index, fontsize=FONT_SIZE)

    ax_heat.set_xticks(np.arange(n_cols))
    pos_labels = []
    for c in mat.columns:
        if c == "Haplogroup":
            pos_labels.append("")
        else:
            try:
                pos = c.split("_")[1]
            except IndexError:
                pos = c
            pos_labels.append(pos)

    tick_labels = []
    for col_name, pos in zip(mat.columns, pos_labels):
        if col_name in sig_variants:
            tick_labels.append(f"{pos}*")
        else:
            tick_labels.append(pos)

    #ax_heat.set_xticklabels(tick_labels, rotation=0, fontsize=FONT_SIZE)
    ax_heat.set_xticks(np.arange(len(tick_labels[1:])) + 1 , tick_labels[1:], rotation=90, fontsize=FONT_SIZE)
    ax_heat.tick_params(axis='x', which='both', pad=2)
    ax_heat.tick_params(axis='y', which='both', pad=2)
    ax_tree.set_axis_off()
    ax_heat.set_xlabel("chr19 position", fontsize=FONT_SIZE)
    ax_heat.spines['top'].set_visible(False)
    ax_heat.spines['right'].set_visible(False)
    ax_heat.spines['bottom'].set_visible(False)
    ax_heat.spines['left'].set_visible(False)

    fig.subplots_adjust(wspace=0.3)

    hap_handles = [Patch(facecolor=hap_colors[hap], edgecolor='k', label=hap) for hap in haplogroups]
    variant_handles = [
        Patch(facecolor="white", edgecolor='k', label="0/0"),
        Patch(facecolor="lightgrey", edgecolor='k', label="0/1"),
        Patch(facecolor="dimgray", edgecolor='k', label="1/1"),
    ]

    hap_legend = fig.legend(
        handles=hap_handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        frameon=False,
        ncol=1,
        fontsize=FONT_SIZE,
        title="Haplogroup",
        title_fontsize=FONT_SIZE
    )

    var_legend = fig.legend(
        handles=variant_handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 0.6),
        frameon=False,
        ncol=1,
        fontsize=FONT_SIZE,
        title="Variant",
        title_fontsize=FONT_SIZE
    )

    fig.add_artist(hap_legend)
    fig.add_artist(var_legend)

    plt.savefig(args.output_pdf, format="pdf", bbox_inches="tight", dpi=600)
    plt.close(fig)

if __name__ == "__main__":
    main()