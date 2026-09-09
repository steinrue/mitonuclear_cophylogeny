configfile: "config.yaml"

wildcard_constraints:
  POP = "|".join(["GWD_YRI", "GWD_PUR", "GWD_CLM", "GWD_CHS", "GWD_JPT", "GWD_IBS", "GWD_TSI", "GWD_GIH", "GWD_STU", "YRI_PUR", "YRI_CLM", "YRI_CHS", "YRI_JPT", "YRI_IBS", "YRI_TSI", "YRI_GIH", "YRI_STU", "PUR_CLM", "PUR_CHS", "PUR_JPT", "PUR_IBS", "PUR_TSI", "PUR_GIH", "PUR_STU", "CLM_CHS", "CLM_JPT", "CLM_IBS", "CLM_TSI", "CLM_GIH", "CLM_STU", "CHS_JPT", "CHS_IBS", "CHS_TSI", "CHS_GIH", "CHS_STU", "JPT_IBS", "JPT_TSI", "JPT_GIH", "JPT_STU", "IBS_TSI", "IBS_GIH", "IBS_STU", "TSI_GIH", "TSI_STU", "GIH_STU", "CEU_GBR", "CEU_KHV", "BEB_MSL", "BEB_ACB", "CHS_CHB", "MSL_ACB", "MSL_GBR", "STU_CHS", "JPT_CHS", "STU_YRI", "CHS_YRI", "GBR", "FIN", "CHS", "PUR", "CDX", "CLM", "IBS", "PEL", "PJL", "KHV", "ACB", "GWD", "ESN", "BEB", "MSL", "STU", "ITU", "CEU", "YRI", "CHB", "JPT", "LWK", "ASW", "MXL", "TSI", "GIH"]),
  CHR = "|".join(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"])

populations_1000G = ["ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI", "CLM", "MXL", "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU"]
paired_chr20 = ["GWD_YRI", "GWD_PUR", "GWD_CLM", "GWD_CHS", "GWD_JPT", "GWD_IBS", "GWD_TSI", "GWD_GIH", "GWD_STU", "YRI_PUR", "YRI_CLM", "YRI_CHS", "YRI_JPT", "YRI_IBS", "YRI_TSI", "YRI_GIH", "YRI_STU", "PUR_CLM", "PUR_CHS", "PUR_JPT", "PUR_IBS", "PUR_TSI", "PUR_GIH", "PUR_STU", "CLM_CHS", "CLM_JPT", "CLM_IBS", "CLM_TSI", "CLM_GIH", "CLM_STU", "CHS_JPT", "CHS_IBS", "CHS_TSI", "CHS_GIH", "CHS_STU", "JPT_IBS", "JPT_TSI", "JPT_GIH", "JPT_STU", "IBS_TSI", "IBS_GIH", "IBS_STU", "TSI_GIH", "TSI_STU", "GIH_STU"]
all_chrs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
chr_20 = ["20"]

rule get_biallelic_vcf:
  input:
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.vcf.gz", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}"),
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.vcf.gz.tbi", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}")
  resources:
    time="2:00:00",
    threads=1,
    mem="10gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}_biallelic.vcf.gz", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}"),
  shell:
    "bcftools view -m2 -M2 -v snps --min-af 0.01 {input[0]} -Oz -o {output} "

rule get_merged_vcf:
  input:
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}_biallelic.vcf.gz", tmp_dir=config['tmp_dir'], POP="{POP}", CHR=all_chrs)
  resources:
    time="1:00:00",
    threads=1,
    mem="10gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}_merged.vcf.gz", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}")
  shell:
    "bcftools concat -Oz -o {output} {input} "

rule get_hardy_stats:
  input:
    expand("{tmp_dir}/{POP}/{POP}_merged.vcf.gz", tmp_dir=config['tmp_dir'], POP="{POP}"),
    expand("{data_dir}/Relate_input_files/GRCh38/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed.gz", data_dir=config['data_dir'], CHR="{CHR}")
  resources:
    time="1:00:00",
    threads=1,
    mem="12gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}.hardy", tmp_dir=config['tmp_dir'], POP="{POP}"),
  shell:
    "plink2 --vcf {input[0]} --extract bed0 {input[1]} "
    "--hardy midp "
    "--out {config[tmp_dir]}/{wildcards.POP}/{wildcards.POP} "
  
rule get_pop_fis:
  input:
    expand("{tmp_dir}/{POP}/{POP}.hardy", tmp_dir=config['tmp_dir'], POP="{POP}")
  resources:
    time="1:00:00",
    threads=1,
    mem="24gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}.fis", tmp_dir=config['tmp_dir'], POP="{POP}"),
  shell:
    "python scripts/fis_jackknife.py {input} > {output}"

rule get_final_fis:
  input:
    expand("{tmp_dir}/{POP}/{POP}.fis", tmp_dir=config['tmp_dir'], POP=populations_1000G),
    expand("{tmp_dir}/{POP}/{POP}.fis", tmp_dir=config['tmp_dir'], POP=paired_chr20)
  resources:
    time="5:00:00",
    threads=1,
    mem="24gb"
  output:
    expand("{results_dir}/FIS/all_populations.fis", results_dir=config['results_dir'])
  shell:
    """
    (echo -e "population\tfis\tlower\tupper"; cat {input}) > {output}
    """
  
rule get_pop_zscore:
  input:
    expand("{results_dir}/all_chrs/{POP}.results", results_dir=config['results_dir'], POP="{POP}")
  resources:
    time="3:00:00",
    threads=1,
    mem="60gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}.zscore", tmp_dir=config['tmp_dir'], POP="{POP}"),
  shell:
    "python scripts/cophylogeny_score_jackknife.py {input} > {output}"

rule get_pop_chr20_zscore:
  input:
    expand("{results_dir}/chr20/{POP}.results", results_dir=config['results_dir'], POP="{POP}")
  resources:
    time="3:00:00",
    threads=1,
    mem="60gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}_chr20.zscore", tmp_dir=config['tmp_dir'], POP="{POP}"),
  shell:
    "python scripts/cophylogeny_score_jackknife.py {input} > {output}"

rule get_final_zscores:
  input:
    expand("{tmp_dir}/{POP}/{POP}.zscore", tmp_dir=config['tmp_dir'], POP=populations_1000G)
  resources:
    time="5:00:00",
    threads=1,
    mem="2gb"
  output:
    expand("{results_dir}/zscores/whole_genome.zscores", results_dir=config['results_dir'])
  shell:
    """
    (echo -e "population\tp_value\tz_score\tlower\tupper"; cat {input}) > {output}
    """

rule get_chr20_zscores:
  input:
    expand("{tmp_dir}/{POP}/{POP}_chr20.zscore", tmp_dir=config['tmp_dir'], POP=paired_chr20)
  resources:
    time="5:00:00",
    threads=1,
    mem="2gb"
  output:
    expand("{results_dir}/zscores/chr20.zscores", results_dir=config['results_dir'])
  shell:
    """
    (echo -e "population\tp_value\tz_score\tlower\tupper"; cat {input}) > {output}
    """
