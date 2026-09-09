import argparse
import pandas as pd
import numpy as np
import os
from scipy.stats import norm
import numba

@numba.njit
def correlateNormals(marginaLoc, marginalScale, corr, numValues, thisRNG, startingValue=None):
    # maybe helpful
    locOne = marginaLoc
    locTwo = marginaLoc
    scaleOne = marginalScale
    scaleTwo = marginalScale

    # get some initial standard normals
    toReturn = thisRNG.normal(loc=0, scale=1, size=numValues)

    # what about first one
    if (startingValue is None):
        # adjust the first one, first scale, then shift
        toReturn[0] = (toReturn[0] * scaleOne) + locOne
    else:
        toReturn[0] = startingValue

    # and then go through to correlate
    for l in np.arange(1, len(toReturn)):
        nextMean = locTwo + (corr * (scaleTwo / scaleOne) * (toReturn[l-1] - locOne))
        nextVar = (1 - corr * corr) * scaleTwo * scaleTwo
        # corrSeq[l] is standard gaussian, first scale, then shift
        toReturn[l] = (toReturn[l] * np.sqrt(nextVar)) + nextMean
    
    return toReturn

def get_pop_ind_thresh(filename):
    basename = os.path.basename(filename)
    base = basename.replace(".txt", "")
    # Expect: {population}_enrichment_{ind}_thresh_{thresh}
    pop_part, rest = base.split('_enrichment_', 1)
    ind_part, thresh_part = rest.split('_thresh_', 1)
    
    population = pop_part
    ind = ind_part
    thresh = float(thresh_part)
    
    return population, ind, thresh

def make_seed_from_strings(population, ind, thresh_str):
    vals = list(population + ind + thresh_str)
    nums = [ord(c) for c in vals]
    return int(np.sum(nums))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--all', required=True, help='All samples pkl')
    parser.add_argument('--output', required=True, help='Output tsv')
    parser.add_argument('--num_replicates', type=int, required=True, help='Number of null genome replicates')
    args = parser.parse_args()

    population, ind, thresh = get_pop_ind_thresh(args.output)
    print(get_pop_ind_thresh(args.output))
    thresh_str = str(thresh)
    seed = make_seed_from_strings(population, ind, thresh_str)
    metaRNG = np.random.default_rng(seed)

    unmasked_df = pd.read_pickle(args.all)

    df = unmasked_df[unmasked_df["in_mask"] == 1].reset_index(drop=True)

    positions = np.array(df['Position'])
    adj_mask = (np.diff(positions // 1000) == 1)
    lp = df[:-1][adj_mask]
    rp = df[1:][adj_mask]

    stepCorr = np.corrcoef(
        norm.isf(10 ** -lp[population]),
        norm.isf(10 ** -rp[population])
    )[0, 1]

    # marginal parameters for correlateNormals
    normalLoc = 0.0
    normalScale = 1.0

    numReplicates = args.num_replicates
    null_vals = []

    all_chroms = unmasked_df["Chromosome"].values
    unique_chroms = unmasked_df["Chromosome"].unique()
    in_mask_array = (unmasked_df["in_mask"].values == 1)

    for rep in np.arange(numReplicates):
        print(f"rep {rep}")
        rep_values = []

        for chr in unique_chroms:
            chr_mask = (all_chroms == chr)
            sequenceLength = int(np.sum(chr_mask))

            # generate correlated normals using the exact function/code you provided
            corrvalues = correlateNormals(normalLoc, normalScale, stepCorr, sequenceLength, metaRNG)
            rep_values.extend(corrvalues)

        rep_values = np.array(rep_values)

        # restrict to in_mask == 1 (same order as df after reset_index)
        rep_vals_masked = rep_values[in_mask_array]

        cutoff = np.quantile(rep_vals_masked, thresh)
        passing = rep_vals_masked >= cutoff  # boolean array, length len(df)

        # subset REAL indicator values using positional mask
        tmp_df = df.loc[passing]

        null_mean = tmp_df[ind].mean()
        null_vals.append(null_mean)

    # observed (real) data: use same positional style
    df_population = df[population].values
    real_cutoff = np.quantile(df_population, thresh)
    real_passing = df_population >= real_cutoff
    real_tmp_df = df.loc[real_passing]

    mean_val = real_tmp_df[ind].mean()
    global_mean = df[ind].mean()

    null_vals_array = np.array(null_vals)
    n = len(null_vals_array)
    p_up = (np.sum(null_vals_array >= mean_val) + 1) / (n + 1)
    p_down = (np.sum(null_vals_array <= mean_val) + 1) / (n + 1)
    pval = min(1.0, 2 * min(p_up, p_down))

    enrichment = mean_val / global_mean

    with open(args.output, 'w') as f_out:
        row = [
            ind,
            population,
            str(thresh),
            str(global_mean),
            str(mean_val),
            str(enrichment),
            str(pval),
            ",".join(str(v) for v in null_vals)
        ]
        f_out.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()