import argparse
import pandas as pd
import numpy as np
import sys
import gzip
from tqdm import tqdm

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
    df = pd.DataFrame(rows, columns=['chr', 'feature', 'start', 'end', 'gene_id', 'gene_name', 'gene_type'])
    df = df[df['gene_type'] == 'protein_coding']
    return df

def main():
    parser = argparse.ArgumentParser(description="Fast annotation of sampled positions with log10_pval from result intervals")
    parser.add_argument('--gff3', required=True, help='Gencode GFF3 file')
    parser.add_argument('--mask', required=True, help='Mask (Bed format)')
    parser.add_argument('--mitocarta', required=True, help='mitocarta file')
    parser.add_argument('--output', required=True, help='Output TSV file')
    args = parser.parse_args()

    output_rows = []
    processed_sample_rows = 0

    mask_df = pd.read_csv(args.mask, sep="\t", header=None, names=['chr', 'start', 'end', 'type'])

    df = parse_gff3_protein_coding_df(args.gff3)
    gene_df = df[df['feature'] == 'gene']
    exon_df = df[df['feature'] == 'exon']
    
    mito_df = pd.read_excel(args.mitocarta, sheet_name=1)
    mito_genes = set()
    for ids in mito_df['EnsemblGeneID_mapping_version_20200130'].dropna().astype(str):
        for gene_id in ids.split('|'):
            mito_genes.add(gene_id.split('.')[0])
    # OXPHOS-specific mito genes
    oxphos_genes = set()
    oxphos_rows = mito_df.loc[mito_df["MitoCarta3.0_MitoPathways"].str.contains("OXPHOS", na=False, regex=False),"EnsemblGeneID_mapping_version_20200130"].dropna().astype(str)
    for ids in oxphos_rows:
        for gene_id in ids.split('|'):
            oxphos_genes.add(gene_id.split('.')[0])
    
    rows = []

    for chrom in CHR_LENGTHS_GRCH38:
        chr_len = CHR_LENGTHS_GRCH38[chrom]
        gene_chr = gene_df[gene_df['chr'] == chrom]
        exon_chr = exon_df[exon_df['chr'] == chrom]
        mask_chr = mask_df[mask_df['chr'] == chrom]
        mask_start = np.array(mask_chr['start'])
        mask_end = np.array(mask_chr['end'])

        sys.stderr.write(f"Processing {chrom} ...\n")
        for pos in tqdm(range(0, chr_len, 1000), desc=f"{chrom}", unit="pos"):
            gene_hits = gene_chr[(gene_chr['start'] <= pos) & (pos < gene_chr['end'])]
            exon_hits = exon_chr[(exon_chr['start'] <= pos) & (pos < exon_chr['end'])]

            mask_idx = np.searchsorted(mask_start, pos, side='right') - 1
            in_mask = 1 if ((pos < mask_end[mask_idx]) & (pos >= mask_start[mask_idx])) else 0
            
            in_gene = 1 if not gene_hits.empty else 0
            gene_ids = ";".join(sorted(gene_hits['gene_id'].unique())) if not gene_hits.empty else ""
            gene_names = ";".join(sorted(gene_hits['gene_name'].unique())) if not gene_hits.empty else ""
            in_exon = 1 if not exon_hits.empty else 0
            exon_ids = ";".join(sorted(exon_hits['gene_id'].unique())) if not exon_hits.empty else ""
            exon_names = ";".join(sorted(exon_hits['gene_name'].unique())) if not exon_hits.empty else ""

            hit_gene_ids = set(gene_hits['gene_id'].unique())
            hit_exon_ids = set(exon_hits['gene_id'].unique())

            mito_gene_overlap = hit_gene_ids & mito_genes
            oxphos_gene_overlap = hit_gene_ids & oxphos_genes

            in_mitocarta_gene = 1 if mito_gene_overlap else 0
            mitocarta_gene_ids = ";".join(sorted(mito_gene_overlap)) if mito_gene_overlap else ""
            mitocarta_gene_names = ";".join(sorted(gene_hits[gene_hits['gene_id'].isin(mito_gene_overlap)]['gene_name'].unique())) if mito_gene_overlap else ""

            in_oxphos_gene = 1 if oxphos_gene_overlap else 0
            oxphos_gene_ids = ";".join(sorted(oxphos_gene_overlap)) if oxphos_gene_overlap else ""
            oxphos_gene_names = ";".join(sorted(gene_hits[gene_hits['gene_id'].isin(oxphos_gene_overlap)]['gene_name'].unique())) if oxphos_gene_overlap else ""

            mito_exon_overlap = hit_exon_ids & mito_genes
            oxphos_exon_overlap = hit_exon_ids & oxphos_genes

            in_mitocarta_exon = 1 if mito_exon_overlap else 0
            mitocarta_exon_ids = ";".join(sorted(mito_exon_overlap)) if mito_exon_overlap else ""
            mitocarta_exon_names = ";".join(sorted(exon_hits[exon_hits['gene_id'].isin(mito_exon_overlap)]['gene_name'].unique())) if mito_exon_overlap else ""

            in_oxphos_exon = 1 if oxphos_exon_overlap else 0
            oxphos_exon_ids = ";".join(sorted(oxphos_exon_overlap)) if oxphos_exon_overlap else ""
            oxphos_exon_names = ";".join(sorted(exon_hits[exon_hits['gene_id'].isin(oxphos_exon_overlap)]['gene_name'].unique())) if oxphos_exon_overlap else ""

            rows.append({
                "Chromosome": chrom,
                "Position": pos,
                "in_mask": in_mask,

                "in_gene": in_gene,
                "gene_ids": gene_ids,
                "gene_names": gene_names,

                "in_exon": in_exon,
                "exon_ids": exon_ids,
                "exon_names": exon_names,

                "in_mitocarta_gene": in_mitocarta_gene,
                "mitocarta_gene_ids": mitocarta_gene_ids,
                "mitocarta_gene_names": mitocarta_gene_names,

                "in_oxphos_gene": in_oxphos_gene,
                "oxphos_gene_ids": oxphos_gene_ids,
                "oxphos_gene_names": oxphos_gene_names,

                "in_mitocarta_exon": in_mitocarta_exon,
                "mitocarta_exon_ids": mitocarta_exon_ids,
                "mitocarta_exon_names": mitocarta_exon_names,

                "in_oxphos_exon": in_oxphos_exon,
                "oxphos_exon_ids": oxphos_exon_ids,
                "oxphos_exon_names": oxphos_exon_names,})

        sys.stderr.write(f"Completed {chrom}.\n")

    result_df = pd.DataFrame(rows, columns=[
        "Chromosome", "Position", "in_mask",

        "in_gene", "gene_ids", "gene_names",
        "in_exon", "exon_ids", "exon_names",

        "in_mitocarta_gene", "mitocarta_gene_ids", "mitocarta_gene_names",
        "in_oxphos_gene", "oxphos_gene_ids", "oxphos_gene_names",
        "in_mitocarta_exon", "mitocarta_exon_ids", "mitocarta_exon_names",
        "in_oxphos_exon", "oxphos_exon_ids", "oxphos_exon_names",])
    result_df.to_pickle(args.output)

if __name__ == "__main__":
    main()