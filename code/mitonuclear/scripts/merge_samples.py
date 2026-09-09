import argparse
import pandas as pd
import os
import sys

def get_pop_name(filename):
    basename = os.path.basename(filename)
    return basename.split('_sampled')[0]

def main():
    parser = argparse.ArgumentParser(description="Merge main loci with population result files by Chromosome/Position")
    parser.add_argument('--gene_info', required=True, help='Main loci annotation file (sampled_loci_phx_mito.tsv)')
    parser.add_argument('--output', required=True, help='Output pickle file')
    parser.add_argument('--popfiles', nargs='+', required=True, help='List of {pop}_sampled.tsv files, space-separated')
    args = parser.parse_args()

    print(f"Reading main file: {args.gene_info}")
    df_main = pd.read_pickle(args.gene_info)

    # Set index for fast joins
    df_main.set_index(['Chromosome','Position'], inplace=True)

    for pfile in args.popfiles:
        pop = get_pop_name(pfile)
        print(f"Merging population file: {pfile} (as column '{pop}')")
        df_pop = pd.read_csv(pfile, sep='\t', dtype={'Chromosome':str, 'Position':int})
        
        if 'log10_pval' not in df_pop.columns:
            sys.exit(f"ERROR: File {pfile} does not contain 'log10_pval' column")
        df_subset = df_pop[['Chromosome', 'Position', 'log10_pval']].copy().rename(columns={'log10_pval': pop})
        df_subset.set_index(['Chromosome','Position'], inplace=True)

        if df_main.index.difference(df_subset.index).size > 0 or df_subset.index.difference(df_main.index).size > 0:
            sys.exit(
                f"ERROR: Chromosome/Position values do not match between main file ({args.main}) and pop file ({pfile}).\n"
                f"Number in main not in pop: {df_main.index.difference(df_subset.index).size}\n"
                f"Number in pop not in main: {df_subset.index.difference(df_main.index).size}"
            )
        # Merge column
        df_main[pop] = df_subset[pop]
        print(f"  ... merged {pop}, shape now: {df_main.shape}")

    print(f"Saving merged dataframe to: {args.output}")
    df_main.reset_index().to_pickle(args.output)
    print("Done.")

if __name__ == "__main__":
    main()