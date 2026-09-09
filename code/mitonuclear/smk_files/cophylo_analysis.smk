configfile: "config.yaml"

wildcard_constraints:
  POP = "|".join(["GWD_YRI", "GWD_PUR", "GWD_CLM", "GWD_CHS", "GWD_JPT", "GWD_IBS", "GWD_TSI", "GWD_GIH", "GWD_STU", "YRI_PUR", "YRI_CLM", "YRI_CHS", "YRI_JPT", "YRI_IBS", "YRI_TSI", "YRI_GIH", "YRI_STU", "PUR_CLM", "PUR_CHS", "PUR_JPT", "PUR_IBS", "PUR_TSI", "PUR_GIH", "PUR_STU", "CLM_CHS", "CLM_JPT", "CLM_IBS", "CLM_TSI", "CLM_GIH", "CLM_STU", "CHS_JPT", "CHS_IBS", "CHS_TSI", "CHS_GIH", "CHS_STU", "JPT_IBS", "JPT_TSI", "JPT_GIH", "JPT_STU", "IBS_TSI", "IBS_GIH", "IBS_STU", "TSI_GIH", "TSI_STU", "GIH_STU", "CEU_GBR", "CEU_KHV", "BEB_MSL", "BEB_ACB", "CHS_CHB", "MSL_ACB", "MSL_GBR", "STU_CHS", "JPT_CHS", "STU_YRI", "CHS_YRI", "GBR", "FIN", "CHS", "PUR", "CDX", "CLM", "IBS", "PEL", "PJL", "KHV", "ACB", "GWD", "ESN", "BEB", "MSL", "STU", "ITU", "CEU", "YRI", "CHB", "JPT", "LWK", "ASW", "MXL", "TSI", "GIH"]),
  CHR = "|".join(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"])

populations_1000G = ["ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI", "CLM", "MXL", "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU"]
paired_chr20 = ["GWD_YRI", "GWD_PUR", "GWD_CLM", "GWD_CHS", "GWD_JPT", "GWD_IBS", "GWD_TSI", "GWD_GIH", "GWD_STU", "YRI_PUR", "YRI_CLM", "YRI_CHS", "YRI_JPT", "YRI_IBS", "YRI_TSI", "YRI_GIH", "YRI_STU", "PUR_CLM", "PUR_CHS", "PUR_JPT", "PUR_IBS", "PUR_TSI", "PUR_GIH", "PUR_STU", "CLM_CHS", "CLM_JPT", "CLM_IBS", "CLM_TSI", "CLM_GIH", "CLM_STU", "CHS_JPT", "CHS_IBS", "CHS_TSI", "CHS_GIH", "CHS_STU", "JPT_IBS", "JPT_TSI", "JPT_GIH", "JPT_STU", "IBS_TSI", "IBS_GIH", "IBS_STU", "TSI_GIH", "TSI_STU", "GIH_STU"]
paired_all_chrs = ["GWD_CHS","GWD_IBS","GWD_PUR","GWD_GIH","CHS_IBS","PUR_CHS","CHS_GIH","PUR_IBS","IBS_GIH","PUR_GIH"]
all_chrs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
chr_20 = ["20"]

rule aggregate_chromosome_20_results:
  input:
    expand("{tmp_dir}/{POP}/results/{POP}_chr{CHR}.results", tmp_dir=config['tmp_dir'], POP="{POP}", CHR=chr_20)
  resources:
    time="0:30:00",
    threads=1,
    mem="3gb"
  output:
    expand("{results_dir}/chr20/{POP}.results", results_dir=config['results_dir'], POP="{POP}")
  shell:
    '''
    bash scripts/aggregate_chromosome_results.sh {input} > {output}
    '''

rule aggregate_chromosome_all_results:
  input:
    expand("{tmp_dir}/{POP}/results/{POP}_chr{CHR}.results", tmp_dir=config['tmp_dir'], POP="{POP}", CHR=all_chrs)
  resources:
    time="0:30:00",
    threads=1,
    mem="3gb"
  output:
    expand("{results_dir}/all_chrs/{POP}.results", results_dir=config['results_dir'], POP="{POP}")
  shell:
    '''
    bash scripts/aggregate_chromosome_results.sh {input} > {output}
    '''

rule get_1kg_mitochondrial_data:
  input:
  resources:
    time="1:00:00",
    threads=1,
    mem="3gb"
  output:
    expand("{data_dir}/1KG_NYGC/1kGP_3202_haplogroup_results.txt", data_dir=config['data_dir'])
  shell:
    "wget -P {config[data_dir]}/1KG_NYGC https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20221104_3202_mitochondrial_pipeline/1kGP_3202_haplogroup_results.txt"

rule get_mitochondria_hsd:
  input:
    expand("{tmp_dir}/{POP}/{POP}_samples.list", tmp_dir=config['tmp_dir'], POP="{POP}"),
    expand("{data_dir}/1KG_NYGC/1kGP_3202_haplogroup_results.txt", data_dir=config['data_dir'])
  resources:
    time="1:00:00",
    threads=1,
    mem="1gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}_MT.hsd", tmp_dir=config['tmp_dir'], POP="{POP}")
  shell:
    "python3 scripts/parse_NYGC_mitochondria.py {input} > {output}"

rule get_mitochondria_msa:
  input:
    expand("{tmp_dir}/{POP}/{POP}_MT.hsd", tmp_dir=config['tmp_dir'], POP="{POP}")
  resources:
    time="1:00:00",
    threads=1,
    mem="10gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}_MT_unclean", tmp_dir=config['tmp_dir'], POP="{POP}"),
    expand("{tmp_dir}/{POP}/{POP}_MT_unclean_MSA.fasta", tmp_dir=config['tmp_dir'], POP="{POP}"),
    expand("{tmp_dir}/{POP}/{POP}_MT_MSA.fasta", tmp_dir=config['tmp_dir'], POP="{POP}")
  run:
    shell("haplogrep3 classify --in {input} --tree phylotree-fu-rcrs@1.2 --output {output[0]} --write-fasta-msa")
    shell("bash scripts/clean_fasta.sh {output[1]} > {output[2]}")

rule get_mitochondria_trees:
  input:
    expand("{tmp_dir}/{POP}/{POP}_MT_MSA.fasta", tmp_dir=config['tmp_dir'], POP="{POP}")
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{trees_dir}/{POP}/{POP}_chrMT.nwk", trees_dir=config['trees_dir'], POP="{POP}"),
    expand("{trees_dir}/{POP}/{POP}_chrMT_shuffled.pkl", trees_dir=config['trees_dir'], POP="{POP}"),
  shell:
    '''
    python3 scripts/get_MT_trees.py -mt_msa {input} -true_tree {output[0]} -shuffled_trees {output[1]} --permutations 1000
    '''

rule aggregate_shuffled_tree_distances:
  input:
    expand("{tmp_dir}/{POP}/tree_distances/{POP}_chr{CHR}_tree_distances.{seeds}", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}", seeds=[str(i) for i in range(0, config['total_parallel'])])
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{tmp_dir}/{POP}/results/{POP}_chr{CHR}.results", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}")
  shell:
    '''
    cat {input} > {output}
    '''

rule get_tree_distances:
  input:
    expand("{tmp_dir}/{POP}/{POP}_samples.list",  tmp_dir=config['tmp_dir'], POP="{POP}"),
    expand("{trees_dir}/{POP}/{POP}_chr{CHR}.trees", trees_dir=config['trees_dir'], POP="{POP}", CHR="{CHR}"),
    expand("{trees_dir}/{POP}/{POP}_chrMT.nwk", trees_dir=config['trees_dir'], POP="{POP}"),
    expand("{trees_dir}/{POP}/{POP}_chrMT_shuffled.pkl", trees_dir=config['trees_dir'], POP="{POP}"),
  resources:
    time="32:00:00",
    threads=1,
    mem="15gb"
  params:
    total_parallel = config["total_parallel"]
  output:
    expand("{tmp_dir}/{POP}/tree_distances/{POP}_chr{CHR}_tree_distances.{SEED}", tmp_dir=config['tmp_dir'], POP="{POP}", SEED ="{SEED}", CHR="{CHR}")
  shell:
    '''
    python3 scripts/tree_distance.py --sample_list_path {input[0]} --tree_sequence_file {input[1]} --mito_tree {input[2]} --shuffled_trees {input[3]} --output_file {output} --parallele_idx {wildcards.SEED} --total_para {params.total_parallel} --permutations 1000 
    '''