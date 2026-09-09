import argparse
import pandas as pd
import numpy as np
import sys
import gzip

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

def main():
    parser = argparse.ArgumentParser(description="Fast annotation of sampled positions with log10_pval from result intervals")
    parser.add_argument('--results', required=True, help='Result intervals TSV')
    parser.add_argument('--output', required=True, help='Output TSV file')
    args = parser.parse_args()

    df_results = pd.read_csv(args.results, sep='\t', dtype={'chr':str})

    df_results['start'] = df_results['start_pos'].apply(lambda x: int(float(x)))

    output_rows = []
    processed_sample_rows = 0

    chroms = sorted(set(df_results['chr']))

    with open(args.output, 'w') as out:
        header = "\t".join(["Chromosome","Position", "log10_pval"])
        print(header, file=out)
        for chrom in CHR_LENGTHS_GRCH38:
            chr_len = CHR_LENGTHS_GRCH38[chrom]
            results_chr = df_results[df_results['chr'] == chrom]
            result_pos = np.array(results_chr['start'])
            result_pval = np.array(results_chr['log10_pval'])

            sys.stderr.write(f"Processing {chrom} ...\n")
            for pos in range(0, chr_len, 1000):
                idx = np.searchsorted(result_pos, pos, side='right') - 1
                pval = result_pval[idx]
                row = [chrom, str(pos), str(pval)]
                print("\t".join(row), file=out)
            sys.stderr.write(f"Completed {chrom}.\n")
    sys.stderr.write(f"All chromosomes processed. Output written to {args.output}\n")

if __name__ == "__main__":
    main()