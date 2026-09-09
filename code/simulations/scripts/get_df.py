import argparse
import pandas as pd
import os
import sys
import numpy as np

def parse_filename(filename):
    base = os.path.basename(filename)
    if base.endswith('.results'):
        name = base[:-len('.results')]
    elif base.endswith('.fis'):
        name = base[:-len('.fis')]
    elif base.endswith('.txt'):
        name = base[:-len('.txt')]
        if name.startswith('allele_freqs_'):
            name = name[len('allele_freqs_'):]
    else:
        return None
    parts = name.split('_')
    params = {}
    if len(parts) % 2 != 0:
        return None
    for key, val in zip(parts[0::2], parts[1::2]):
        try:
            if val.isdigit():
                val_cast = int(val)
            else:
                val_cast = float(val)
        except ValueError:
            val_cast = val
        params[key] = val_cast
    return params

def main():
    parser = argparse.ArgumentParser(description="Combine .results files, optionally with allele frequency files. Merges and adds filename parameters as columns.")
    parser.add_argument("-output_file", required=True, help="Path to the output pickle file for the combined DataFrame")
    parser.add_argument("-input_files", required=True, nargs="+",help="Space-separated list of .results and optionally allele frequency .txt files")
    parser.add_argument("--allele_freq_data", action="store_true", help="If set, will parse allele_freqs_*.txt files and merge with results",default=False)
    args = parser.parse_args()

    results_dfs = []
    allele_freq_dfs = []

    for f in args.input_files:
        if f.endswith('.results'):
            params = parse_filename(f)
            if not params:
                continue
            df = pd.read_csv(f, sep="\t", header=0)
            for k, v in params.items():
                df[k] = v
            results_dfs.append(df)
        elif f.endswith('.fis'):
            params = parse_filename(f)
            if not params:
                continue
            df = pd.read_csv(f, sep="\t", header=0)
            for k, v in params.items():
                df[k] = v
            results_dfs.append(df)
        elif args.allele_freq_data and f.endswith('.txt') and 'allele_freqs_' in os.path.basename(f):
            params = parse_filename(f)
            if not params:
                continue
            freq_df = pd.read_csv(f, sep=",", header=0)
            freq_df['tick'] = freq_df['tick'].astype(int)
            freq_df['ST'] = freq_df['tick'] - 150000
            freq_df = freq_df.drop(columns=['tick'])
            for k, v in params.items():
                freq_df[k] = v
            allele_freq_dfs.append(freq_df)

    if not results_dfs:
        print("Error: No valid result files provided.", file=sys.stderr)
        sys.exit(1)
    
    results_all = pd.concat(results_dfs, ignore_index=True)

    if args.allele_freq_data:
        if allele_freq_dfs:
            allele_freq_all = pd.concat(allele_freq_dfs, ignore_index=True)
            allele_freq_all = allele_freq_all[allele_freq_all["ST"] == 1]
            if 'PROP' in params.keys():
                merge_on = ['PM', 'PR','PROP', 'S1', 'S2', 'seed']
            else:
                merge_on = ['PM', 'PR', 'S1', 'S2', 'seed']

            merged = pd.merge(
                results_all,
                allele_freq_all[merge_on + ['host_A', 'host_A_loc', 'symb_A']],
                on=merge_on,
                how="left"
            )
        else:
            print("Warning: No valid allele frequency files provided.", file=sys.stderr)
            merged = results_all.copy()
            merged['host_A'] = np.nan
            merged['host_A_loc'] = np.nan
            merged['symb_A'] = np.nan
    else:
        merged = results_all

    merged.to_pickle(args.output_file)
    print(f"Combined DataFrame saved as pickle to {args.output_file}")

if __name__ == "__main__":
    main()