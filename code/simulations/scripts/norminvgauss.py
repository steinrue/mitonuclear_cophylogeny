#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
from scipy.stats import norminvgauss, beta
import pickle

def main():
    parser = argparse.ArgumentParser(description="Compute expected p-values, fitted NIG p-values, and beta-based bounds from a results file.")
    parser.add_argument("--input","-i",required=True,help="Path to the input .results file (tab-separated, with a 'td' column).")
    parser.add_argument("--output","-o",required=True,help="Path to the output pickle file (DataFrame will be saved here).")
    parser.add_argument("--n_rows","-n",type=int,default=None,help="Number of rows from the input file to process (default: all rows).")
    parser.add_argument("--subsample_size","-s",type=int,default=1000,help="Subsample size per row when fitting NIG (default: 1000).")
    args = parser.parse_args()
    path = args.input
    out_path = args.output
    n_rows = args.n_rows
    subsample_size = args.subsample_size
    df = pd.read_csv(path, sep="\t")
    if n_rows is None:
        n_rows = df.shape[0]
    else:
        n_rows = min(n_rows, df.shape[0])
    vals0 = np.array([float(x) for x in df.loc[0,"td"].split(";")])
    n = len(vals0)
    i = np.arange(1, n + 1)
    m = n + 1
    expected_p = i / m
    lower = -np.log10(beta.ppf(0.025, i, m - i))
    upper = -np.log10(beta.ppf(0.975, i, m - i))
    out = pd.DataFrame({"expected_p": expected_p, "lower": lower, "upper": upper})
    for row in range(n_rows):
        v = np.array([float(x) for x in df.loc[row,"td"].split(";")])
        samp = np.random.choice(v, size=subsample_size, replace=False)
        nig_params = norminvgauss.fit(samp)
        p_nig = norminvgauss.cdf(v, *nig_params)
        out[str(row)] = p_nig
    with open(out_path,"wb") as f:
        pickle.dump(out,f)

if __name__ == "__main__":
    main()