import pandas as pd
import numpy as np
import argparse
import os

def load_hardy_file(filepath):
    """Load .hardy file into a DataFrame with clean headers."""
    with open(filepath) as f:
        header = f.readline().lstrip("#").strip().split()
    df = pd.read_csv(filepath, delim_whitespace=True, skiprows=1, names=header)
    return df

# according to Busing et al. (1999): Delete-m Jackknife for Unequal m
def jackknifeBootstrap (wholeTheta, numDatapoints, blockThetas, blockLengths):
    assert (len(blockThetas) == len(blockLengths))
    assert (blockLengths.sum() == numDatapoints)
    numBlocks = len(blockThetas)
    fracBlockSizes = blockLengths / numDatapoints
    # bias correct estimate
    unbiasedTheta = numBlocks * wholeTheta - ((1 - fracBlockSizes)*blockThetas).sum()
    # variance
    hs = 1/fracBlockSizes
    pseudoValues = (hs * wholeTheta) - ((hs - 1)*blockThetas)
    bootstrapVariance = ((pseudoValues - unbiasedTheta)**2/(hs-1)).sum()/numBlocks
    bootStrapSd = np.sqrt(bootstrapVariance)
    return (unbiasedTheta, unbiasedTheta - 1.96*bootStrapSd, unbiasedTheta + 1.96*bootStrapSd)

def ratio_of_averages(df, num_samples, het_obs_col="O(HET_A1)", het_exp_col="E(HET_A1)"):
    sum_obs = df[het_obs_col].sum()
    sum_exp = (2*num_samples/(2*num_samples - 1)) * df[het_exp_col].sum()
    return 1 - (sum_obs / sum_exp)

def jackknife_by_chromosome(df, num_samples, chrom_col="CHROM", het_obs_col="O(HET_A1)", het_exp_col="E(HET_A1)", maf=0.01):
    df = df[df[het_exp_col] > (2*(maf)*(1-maf))].copy()
    roa = ratio_of_averages(df, num_samples)
    numDatapoints = df.shape[0]
    chromosomes = df[chrom_col].unique()
    blockThetas = []
    blockLengths = []
    for chrom in chromosomes:
        df_jack = df[df[chrom_col] != chrom]
        fis = ratio_of_averages(df_jack, num_samples)
        blockThetas.append(fis)
        blockLengths.append(df[df[chrom_col] == chrom].shape[0])
    blockThetas = np.array(blockThetas)
    blockLengths = np.array(blockLengths)
    return jackknifeBootstrap (roa, numDatapoints, blockThetas, blockLengths)

def extract_population_name(filepath):
    return os.path.basename(os.path.dirname(filepath))

def main():
    parser = argparse.ArgumentParser(description="Compute FIS Ratio of Averages + jackknife SE for multiple populations")
    parser.add_argument("hardy_files", nargs="+", help="Paths to .hardy files, one per population")
    args = parser.parse_args()

    for filepath in args.hardy_files:
        pop = extract_population_name(filepath)
        df = load_hardy_file(filepath)

        if "CHROM" not in df.columns:
            df = df.rename(columns={"#CHROM": "CHROM"})
            
        cols = ["HOM_A1_CT", "HET_A1_CT", "TWO_AX_CT"]
        n1 = df.loc[0, cols].sum()
        n2 = df.loc[1, cols].sum()
        assert n1 == n2
        num_samples = int(n1)
        
        theta, lower, upper = jackknife_by_chromosome(df, num_samples)

        print(f"{pop}\t{theta}\t{lower}\t{upper}")

if __name__ == "__main__":
    main()
