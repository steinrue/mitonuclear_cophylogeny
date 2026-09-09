configfile: "config.yaml"

wildcard_constraints:
  POP = "|".join(["GWD_YRI", "GWD_PUR", "GWD_CLM", "GWD_CHS", "GWD_JPT", "GWD_IBS", "GWD_TSI", "GWD_GIH", "GWD_STU", "YRI_PUR", "YRI_CLM", "YRI_CHS", "YRI_JPT", "YRI_IBS", "YRI_TSI", "YRI_GIH", "YRI_STU", "PUR_CLM", "PUR_CHS", "PUR_JPT", "PUR_IBS", "PUR_TSI", "PUR_GIH", "PUR_STU", "CLM_CHS", "CLM_JPT", "CLM_IBS", "CLM_TSI", "CLM_GIH", "CLM_STU", "CHS_JPT", "CHS_IBS", "CHS_TSI", "CHS_GIH", "CHS_STU", "JPT_IBS", "JPT_TSI", "JPT_GIH", "JPT_STU", "IBS_TSI", "IBS_GIH", "IBS_STU", "TSI_GIH", "TSI_STU", "GIH_STU", "CEU_GBR", "CEU_KHV", "BEB_MSL", "BEB_ACB", "CHS_CHB", "MSL_ACB", "MSL_GBR", "STU_CHS", "JPT_CHS", "STU_YRI", "CHS_YRI", "GBR", "FIN", "CHS", "PUR", "CDX", "CLM", "IBS", "PEL", "PJL", "KHV", "ACB", "GWD", "ESN", "BEB", "MSL", "STU", "ITU", "CEU", "YRI", "CHB", "JPT", "LWK", "ASW", "MXL", "TSI", "GIH"]),
  CHR = "|".join(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]),
  IND = "|".join(["in_exon", "in_mitocarta_exon", "in_oxphos_exon", "in_gene", "in_mitocarta_gene", "in_oxphos_gene"]),
  THRESH = "|".join(["0.99"]),

populations_1000G = ["ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI", "CLM", "MXL", "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU"]
paired_all_chrs = ["GWD_CHS","GWD_IBS","GWD_PUR","GWD_GIH","CHS_IBS","PUR_CHS","CHS_GIH","PUR_IBS","IBS_GIH","PUR_GIH"]
all_chrs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]

rule get_gencode_data:
  input:
  resources:
    time="1:00:00",
    threads=1,
    mem="1gb"
  output:
    expand("{data_dir}/GenCode/gencode.v49.annotation.gff3.gz", data_dir=config['data_dir'])
  shell:
    "wget -P {config[data_dir]}/GenCode https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gff3.gz"

rule get_mitocarta_data:
  input:
  resources:
    time="1:00:00",
    threads=1,
    mem="1gb"
  output:
    expand("{data_dir}/MitoCarta/Human.MitoCarta3.0.xls", data_dir=config['data_dir'])
  shell:
    "wget -P {config[data_dir]}/MitoCarta https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls"

rule merge_gene_information:
  input:
    expand("{data_dir}/MitoCarta/Human.MitoCarta3.0.xls", data_dir=config['data_dir']),
    expand("{data_dir}/GenCode/gencode.v49.annotation.gff3.gz", data_dir=config['data_dir']),
    expand("{data_dir}/Relate_input_files/GRCh38/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed.gz", data_dir=config['data_dir'])
  resources:
    time="2:00:00",
    threads=1,
    mem="10gb"
  output:
    expand("{results_dir}/single_locus/sample_information.pkl", POP=populations_1000G, results_dir=config['results_dir']),
  shell:
    "python scripts/gene_sample.py --mitocarta {input[0]} --gff3 {input[1]} --mask {input[2]} --output {output}"

rule merge:
  input:
    gene_info = expand("{results_dir}/single_locus/sample_information.pkl",results_dir=config['results_dir']),
    popfiles = expand("{results_dir}/single_locus/{POP}_sampled.tsv", POP=populations_1000G, results_dir=config['results_dir'])
  resources:
    time="1:10:00",
    threads=1,
    mem="20gb"
  output:
    "{results_dir}/single_locus/all_samples.pkl"
  shell:
    r"""
    python scripts/merge_samples.py \
      --gene_info {input.gene_info[0]} \
      --output {output} \
      --popfiles {input.popfiles}
    """

rule reduce_results:
  input:
    expand("{results_dir}/all_chrs/{POP}.results",  POP="{POP}", results_dir=config['results_dir']),
  resources:
    time="0:10:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/single_locus/{POP}.pvals", POP="{POP}", results_dir=config['results_dir'])
  run:
    shell(
      "awk 'BEGIN{{OFS=\"\t\"}} {{NF--; print}}' {input} > {output}"
    )

rule sample_results:
  input:
    expand("{results_dir}/single_locus/{POP}.pvals", POP="{POP}", results_dir=config['results_dir'])
  resources:
    time="0:10:00",
    threads=1,
    mem="10gb"
  output:
    expand("{results_dir}/single_locus/{POP}_sampled.tsv", POP="{POP}", results_dir=config['results_dir']),
  shell:
    '''
    python scripts/sample_results.py --results {input} --output {output}
    '''

rule pop_enrichment:
  input:
    expand("{results_dir}/single_locus/all_samples.pkl", results_dir=config['results_dir'])
  resources:
    time="5:00:00",
    threads=1,
    mem="12gb"
  output:
    expand("{tmp_dir}/single_locus/{POP}_enrichment_{IND}_thresh_{THRESH}.txt", POP="{POP}", IND="{IND}", THRESH="{THRESH}", tmp_dir=config['tmp_dir']),
  shell:
    '''
    python scripts/enrichment.py --all {input} --num_replicates 10000 --output {output}
    '''

rule enrichment:
  input:
    expand("{tmp_dir}/single_locus/{POP}_enrichment_{IND}_thresh_{THRESH}.txt", POP=populations_1000G, IND=["in_exon", "in_mitocarta_exon", "in_oxphos_exon", "in_gene", "in_mitocarta_gene", "in_oxphos_gene"], THRESH="0.99", tmp_dir=config['tmp_dir']),
  resources:
    time="3:00:00",
    threads=1,
    mem="12gb"
  output:
    expand("{results_dir}/single_locus/enrichment_df.pkl", results_dir=config['results_dir']),
  shell:
      r'''
      cat {input} | python -c "import sys, pandas as pd; df = pd.read_csv(sys.stdin, sep='\t', header=None, names=['ind','population','thresh','global_mean','mean','enrichment','pval','null_vals']); df.to_pickle('{output}')"
      '''