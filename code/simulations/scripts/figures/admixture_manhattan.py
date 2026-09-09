import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import met_brewer


ST_ORDER = [5, 10, 25, 100]
SUBSAMPLE_STEP = 1000
WINDOW = 5e5
FONT_SIZE = 7
FIG_WIDTH_IN = 160 / 25.4
FIG_HEIGHT_IN = 40 / 25.4


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--admixture-pkl",
        required=True,
        help=("Path to admixture.pkl."),
    )
    parser.add_argument(
        "--output-pdf",
        required=True,
        help=("Path to output PDF file."),
    )
    args = parser.parse_args()

    df = pd.read_pickle(args.admixture_pkl)
    df = df[(df["S1"] == 0.0) & (df["S2"] == 1)]
    df = df[df["PROP"] == 0.5]

    colors = met_brewer.met_brew(name="Redon", n=50, brew_type="continuous")
    param_groups = df.groupby(["S1", "S2"])
    group_keys = list(param_groups.groups.keys())
    n_param_rows = len(group_keys)

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
        n_param_rows,
        len(ST_ORDER),
        figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
        sharex=True,
        dpi=600,
        constrained_layout=True,
    )

    if n_param_rows == 1:
        axs = [axs]

    for row_idx, (s1, s2) in enumerate(group_keys):
        group_df = param_groups.get_group((s1, s2))
        seeds = sorted(group_df["seed"].unique())
        color_dict = {}
        for i in range(50):
            if i < len(seeds):
                color_dict[seeds[i]] = colors[i]
        for col_idx, st in enumerate(ST_ORDER):
            ax = axs[row_idx][col_idx] if n_param_rows > 1 else axs[0][col_idx]
            focal_arrays = []
            x_positions = np.arange(-WINDOW, WINDOW + 1, SUBSAMPLE_STEP)
            for seed in seeds:
                query = (group_df["seed"] == seed) & (group_df["ST"] == st)
                if not np.any(query):
                    continue
                row = group_df[query].iloc[0]
                pos = np.array([float(x) for x in str(row["start_vals"]).split(";")])
                vals = np.array([float(x) for x in str(row["p_vals"]).split(";")])
                host_loc = float(row["host_A_loc"])
                sample_positions = host_loc + x_positions
                focal_vals = []
                for sp in sample_positions:
                    indices = np.where(pos <= sp)[0]
                    if len(indices) > 0:
                        idx = indices[-1]
                        focal_vals.append(vals[idx])
                    else:
                        focal_vals.append(np.nan)
                focal_arrays.append(focal_vals)
                ax.scatter(
                    x_positions,
                    focal_vals,
                    color=color_dict[seed],
                    alpha=0.1,
                    s=0.2,
                    rasterized=True
                )

            if len(focal_arrays) > 0:
                median_vals = np.nanmedian(focal_arrays, axis=0)
                ax.plot(
                    x_positions,
                    median_vals,
                    color="k",
                    alpha=0.85,
                    linewidth=1,
                    rasterized=True
                )

            ax.set_title(f"Generation {st}", fontsize=FONT_SIZE)
            ax.axhline(0, color="grey", linestyle="--", linewidth=1)
            ax.scatter(0, 0, color="white", alpha=0)

            if col_idx == 0:
                ax.set_ylabel(r"$-\log_{10}(p)$", fontsize=FONT_SIZE)
            else:
                ax.set_xlabel("", fontsize=FONT_SIZE)

            ax.tick_params(axis="both", labelsize=FONT_SIZE)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
    fig.supxlabel('Base pairs from focal loci', fontsize=FONT_SIZE)
    fig.savefig(args.output_pdf, format="pdf", dpi = 600)
    plt.close(fig)


if __name__ == "__main__":
    main()