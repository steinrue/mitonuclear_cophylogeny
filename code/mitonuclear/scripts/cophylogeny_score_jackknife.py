import pandas as pd
import numpy as np
import argparse
import os

def analyze_td_distribution(df):
    td_nans = df['td'].isna().sum()
    td_total = len(df)
    td_nan_frac = td_nans / td_total
    
    if td_nan_frac > 0.05:
        true_mean = np.nan
        true_median = np.nan
    else:
        true_mean = df['td'].mean(skipna=True)
        true_median = df['td'].median(skipna=True)
    
    pval_nans = df['log10_pval'].isna().sum()
    pval_nan_frac = pval_nans / td_total
    if pval_nan_frac > 0.05:
        log10_pval_mean = np.nan
        log10_pval_median = np.nan
    else:
        log10_pval_mean = df['log10_pval'].mean(skipna=True)
        log10_pval_median = df['log10_pval'].median(skipna=True)

    shuffled_matrix = df['shuffled_td_values'].apply(lambda x: np.array([float(yi) if yi != 'nan' else np.nan for yi in x.strip().split(';')]))
    assert all(len(vals) == 1000 for vals in shuffled_matrix), f"File {filepath} has rows with != 1000 shuffled values"
    
    shuffled_array = np.stack(shuffled_matrix.values)
    nan_fraction_per_col = np.mean(np.isnan(shuffled_array), axis=0)
    valid_cols = nan_fraction_per_col <= 0.05

    percent_invalid = 100 * (np.sum(~valid_cols) / 1000)
    if percent_invalid > 5: print(f"{percent_invalid:.2f}% of shuffled columns are invalid (>5% NANs) and excluded.")
    valid_shuffled = shuffled_array[:, valid_cols]

    if valid_shuffled.shape[1] < 990:
        null_mean = np.nan
        null_median = np.nan
        null_std = np.nan
        p_value = np.nan
    else:
        shuffled_means = np.nanmean(valid_shuffled, axis=0)
        shuffled_median = np.nanmedian(shuffled_means)
        null_mean = np.nanmean(shuffled_means)
        null_median = shuffled_median
        null_std = np.nanstd(shuffled_means, ddof=1)
        if np.isnan(true_mean):
            p_value = np.nan
        else:
            p_value = (np.sum(shuffled_means <= true_mean) + 1) / (len(shuffled_means) + 1)
    
    if any(np.isnan([null_mean, null_std, true_mean])):
        z_score = np.nan
    else:
        z_score = (true_mean - null_mean) / null_std

    return p_value, z_score
    
# according to Busing et al. (1999): Delete-m Jackknife for Unequal m
def jackknifeBootstrap (pval, wholeTheta, numDatapoints, blockThetas, blockLengths):
    assert (len(blockThetas) == len(blockLengths))
    assert (blockLengths.sum() == numDatapoints)
    numBlocks = len(blockThetas)
    fracBlockSizes = blockLengths / numDatapoints
    unbiasedTheta = numBlocks * wholeTheta - ((1 - fracBlockSizes)*blockThetas).sum()
    hs = 1/fracBlockSizes
    pseudoValues = (hs * wholeTheta) - ((hs - 1)*blockThetas)
    bootstrapVariance = ((pseudoValues - unbiasedTheta)**2/(hs-1)).sum()/numBlocks
    bootStrapSd = np.sqrt(bootstrapVariance)
    return (pval, unbiasedTheta, unbiasedTheta - 1.96*bootStrapSd, unbiasedTheta + 1.96*bootStrapSd)

def jackknife_by_chromosome(df, chrom_col="chr"):
    pval, true_score = analyze_td_distribution(df)
    numDatapoints = df.shape[0]
    chromosomes = df[chrom_col].unique()
    blockThetas = []
    blockLengths = []
    for chrom in chromosomes:
        df_jack = df[df[chrom_col] != chrom]
        jack_pval, z_score = analyze_td_distribution(df_jack)
        blockThetas.append(z_score)
        blockLengths.append(df[df[chrom_col] == chrom].shape[0])
    blockThetas = np.array(blockThetas)
    blockLengths = np.array(blockLengths)
    return jackknifeBootstrap (pval, true_score, numDatapoints, blockThetas, blockLengths)

def jackknife_by_chunks(df, n_chunks=20): 
    pval, true_score = analyze_td_distribution(df) 
    numDatapoints = df.shape[0]
    indices = np.arange(numDatapoints)
    chunks = np.array_split(indices, n_chunks)
    blockThetas = []
    blockLengths = []
    for chunk in chunks:
        mask = np.ones(numDatapoints, dtype=bool)
        mask[chunk] = False
        df_jack = df.iloc[mask]
        jack_pval, z_score = analyze_td_distribution(df_jack)
        blockThetas.append(z_score)
        blockLengths.append(len(chunk))
    blockThetas = np.array(blockThetas)
    blockLengths = np.array(blockLengths)
    return jackknifeBootstrap(pval, true_score, numDatapoints, blockThetas, blockLengths)

def main():
    parser = argparse.ArgumentParser(description="Compute FIS Ratio of Averages + jackknife SE for multiple populations")
    parser.add_argument("result_files", nargs="+", help="Paths to .results files, one per population")
    args = parser.parse_args()

    for filepath in args.result_files:
        pop = os.path.basename(filepath).replace('.results', '')
        df = pd.read_csv(filepath, sep='\t')
        if len(df["chr"].unique()) == 1:
            pval, theta, lower, upper = jackknife_by_chunks(df)
        else:
            pval, theta, lower, upper = jackknife_by_chromosome(df)

        print(f"{pop}\t{pval}\t{theta}\t{lower}\t{upper}")

if __name__ == "__main__":
    main()
