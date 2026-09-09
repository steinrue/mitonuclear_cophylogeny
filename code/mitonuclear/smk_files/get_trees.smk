configfile: "config.yaml"

wildcard_constraints:
  POP = "|".join(["GWD_YRI", "GWD_PUR", "GWD_CLM", "GWD_CHS", "GWD_JPT", "GWD_IBS", "GWD_TSI", "GWD_GIH", "GWD_STU", "YRI_PUR", "YRI_CLM", "YRI_CHS", "YRI_JPT", "YRI_IBS", "YRI_TSI", "YRI_GIH", "YRI_STU", "PUR_CLM", "PUR_CHS", "PUR_JPT", "PUR_IBS", "PUR_TSI", "PUR_GIH", "PUR_STU", "CLM_CHS", "CLM_JPT", "CLM_IBS", "CLM_TSI", "CLM_GIH", "CLM_STU", "CHS_JPT", "CHS_IBS", "CHS_TSI", "CHS_GIH", "CHS_STU", "JPT_IBS", "JPT_TSI", "JPT_GIH", "JPT_STU", "IBS_TSI", "IBS_GIH", "IBS_STU", "TSI_GIH", "TSI_STU", "GIH_STU", "CEU_GBR", "CEU_KHV", "BEB_MSL", "BEB_ACB", "CHS_CHB", "MSL_ACB", "MSL_GBR", "STU_CHS", "JPT_CHS", "STU_YRI", "CHS_YRI", "GBR", "FIN", "CHS", "PUR", "CDX", "CLM", "IBS", "PEL", "PJL", "KHV", "ACB", "GWD", "ESN", "BEB", "MSL", "STU", "ITU", "CEU", "YRI", "CHB", "JPT", "LWK", "ASW", "MXL", "TSI", "GIH"]),
  CHR = "|".join(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"])

populations_1000G = ["ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI", "CLM", "MXL", "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU"]
paired_chr20 = ["GWD_YRI", "GWD_PUR", "GWD_CLM", "GWD_CHS", "GWD_JPT", "GWD_IBS", "GWD_TSI", "GWD_GIH", "GWD_STU", "YRI_PUR", "YRI_CLM", "YRI_CHS", "YRI_JPT", "YRI_IBS", "YRI_TSI", "YRI_GIH", "YRI_STU", "PUR_CLM", "PUR_CHS", "PUR_JPT", "PUR_IBS", "PUR_TSI", "PUR_GIH", "PUR_STU", "CLM_CHS", "CLM_JPT", "CLM_IBS", "CLM_TSI", "CLM_GIH", "CLM_STU", "CHS_JPT", "CHS_IBS", "CHS_TSI", "CHS_GIH", "CHS_STU", "JPT_IBS", "JPT_TSI", "JPT_GIH", "JPT_STU", "IBS_TSI", "IBS_GIH", "IBS_STU", "TSI_GIH", "TSI_STU", "GIH_STU"]
#paired_all_chrs = ["GWD_CHS","GWD_IBS","GWD_PUR","GWD_GIH","CHS_IBS","PUR_CHS","CHS_GIH","PUR_IBS","IBS_GIH","PUR_GIH"]
all_chrs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
chr_20 = ["20"]

rule get_input_data:
  input:
  resources:
    time="1:00:00",
    threads=1,
    mem="1gb"
  output:
    expand("{data_dir}/Relate_input_files/coal_rates/1000G_auto.coal", data_dir=config['data_dir']),
    expand("{data_dir}/Relate_input_files/GRCh38/human_ancestor_GRCh38/homo_sapiens_ancestor_{CHR}.fa.gz", data_dir=config['data_dir'], CHR=all_chrs),
    expand("{data_dir}/Relate_input_files/GRCh38/20160622_genome_mask_GRCh38/StrictMask/20160622.chr{CHR}.mask.fasta.gz", data_dir=config['data_dir'], CHR=all_chrs),
    expand("{data_dir}/Relate_input_files/GRCh38/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed.gz", data_dir=config['data_dir']),
    expand("{data_dir}/Relate_input_files/GRCh38/recomb_map/genetic_map_chr{CHR}.txt", data_dir=config['data_dir'], CHR=all_chrs),
  run:
    shell("wget -O {config[data_dir]}/Relate_input_files.tgz https://zenodo.org/records/15801307/files/Relate_input_files.tgz")
    shell("tar -xzf {config[data_dir]}/Relate_input_files.tgz -C {config[data_dir]}")

rule get_1kg_metadata:
  input:
  resources:
    time="1:00:00",
    threads=1,
    mem="1gb"
  output:
    expand("{data_dir}/1KG_NYGC/1000G_2504_high_coverage.sequence.index", data_dir=config['data_dir'])
  shell:
    "wget -P {config[data_dir]}/1KG_NYGC https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index"

rule get_superpop_metadata:
  input:
  resources:
    time="1:00:00",
    threads=1,
    mem="1gb"
  output:
    expand("{data_dir}/1KG_NYGC/20130606_g1k_3202_samples_ped_population.txt", data_dir=config['data_dir'])
  shell:
    "wget -P {config[data_dir]}/1KG_NYGC https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"

rule clean_superpop_metadata:
  input:
    expand("{data_dir}/1KG_NYGC/20130606_g1k_3202_samples_ped_population.txt", data_dir=config['data_dir'])
  resources:
    time="1:00:00",
    threads=1,
    mem="1gb"
  output:
    expand("{data_dir}/1KG_NYGC/1KG.poplabels", data_dir=config['data_dir'])
  shell:
    "python scripts/get_1KG_poplabels.py {input} > {output}"

rule get_1KG_NYGC_data:
  input:
  resources:
    time="2:00:00",
    threads=1,
    mem="1gb"
  output:
    expand("{data_dir}/1KG_NYGC/1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz", data_dir=config['data_dir'], CHR="{CHR}"),
    expand("{data_dir}/1KG_NYGC/1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi", data_dir=config['data_dir'], CHR="{CHR}")
  shell:
    "wget -P {config[data_dir]}/1KG_NYGC https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr{wildcards.CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr{wildcards.CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

rule get_pop_samples:
  input:
    expand("{data_dir}/1KG_NYGC/1000G_2504_high_coverage.sequence.index", data_dir=config['data_dir'])
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}_samples.list",  tmp_dir=config['tmp_dir'], POP="{POP}")
  shell:
    "python scripts/get_sample_list.py {wildcards.POP} {input} > {output}"

rule parse_vcf_file:
  input:
    expand("{tmp_dir}/{POP}/{POP}_samples.list", tmp_dir=config['tmp_dir'], POP="{POP}"),
    expand("{data_dir}/1KG_NYGC/1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz", data_dir=config['data_dir'], CHR="{CHR}"),
    expand("{data_dir}/1KG_NYGC/1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi", data_dir=config['data_dir'], CHR="{CHR}")
  resources:
    time="2:00:00",
    threads=1,
    mem="2gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.vcf", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}")
  run:
    shell("python scripts/parse_vcf_populations.py {input[0]} {input[1]} {output}")

rule zip_and_index_vcf:
  input:
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.vcf", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}")
  resources:
    time="1:00:00",
    threads=1,
    mem="3gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.vcf.gz", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}"),
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.vcf.gz.tbi", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}")
  run:
    shell("bgzip {input}")
    shell("bcftools index -t {output[0]}")

rule get_single_poplabel:
  input:
    expand("{tmp_dir}/{POP}/{POP}_samples.list", tmp_dir=config['tmp_dir'], POP="{POP}"),
    expand("{data_dir}/1KG_NYGC/1KG.poplabels", data_dir=config['data_dir'])
  resources:
    time="1:00:00",
    threads=1,
    mem="1gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}.poplabels", tmp_dir=config['tmp_dir'], POP="{POP}")
  shell:
    "python scripts/get_single_poplabel.py {input} > {output}"

rule vcf_to_hapsamp:
  input:
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.vcf.gz", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}"),
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.vcf.gz.tbi", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}")
  resources:
    time="2:00:00",
    threads=1,
    mem="6gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}_unprepared.haps", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}"),
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}_unprepared.sample", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}")
  shell:
    "{config[path_to_relate]}/bin/RelateFileFormats --mode ConvertFromVcf "
    "--haps {output[0]} "
    "--sample {output[1]} "
    "-i {config[tmp_dir]}/{wildcards.POP}/{wildcards.POP}_chr{wildcards.CHR}/{wildcards.POP}_chr{wildcards.CHR}"

rule prepare_hapsamp:
  input:
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}_unprepared.haps", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}"),
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}_unprepared.sample", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}"),
    expand("{data_dir}/Relate_input_files/GRCh38/human_ancestor_GRCh38/homo_sapiens_ancestor_{CHR}.fa.gz", data_dir=config['data_dir'], CHR="{CHR}"),
    expand("{data_dir}/Relate_input_files/GRCh38/20160622_genome_mask_GRCh38/StrictMask/20160622.chr{CHR}.mask.fasta.gz", data_dir=config['data_dir'], CHR="{CHR}")
  resources:
    time="1:00:00",
    threads=1,
    mem="6gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.haps.gz", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}"),
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.sample.gz", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}"),
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.dist", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}")
  run:
    shell("{config[path_to_relate]}/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps {input[0]} --sample {input[1]} --ancestor {input[2]} --mask {input[3]} -o {config[tmp_dir]}/{wildcards.POP}/{wildcards.POP}_chr{wildcards.CHR}/{wildcards.POP}_chr{wildcards.CHR}")
    shell("gunzip {config[tmp_dir]}/{wildcards.POP}/{wildcards.POP}_chr{wildcards.CHR}/{wildcards.POP}_chr{wildcards.CHR}.dist.gz")

rule run_relate:
  input:
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.haps.gz", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}"),
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.sample.gz", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}"),
    expand("{data_dir}/Relate_input_files/GRCh38/recomb_map/genetic_map_chr{CHR}.txt", data_dir=config['data_dir'], CHR="{CHR}"),
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.dist", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}"),
    expand("{data_dir}/Relate_input_files/coal_rates/1000G_auto.coal", data_dir=config['data_dir']),
  resources:
    time="12:00:00",
    threads=1,
    mem="12gb"
  output:
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.anc", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}"),
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.mut", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}")
  run:
    shell("{config[path_to_relate]}/bin/Relate --mode All -m 1.25e-8 -N 30000 --haps {input[0]} --sample {input[1]} --map {input[2]} --dist {input[3]} --coal {input[4]} --seed 77 --memory 7 -o {wildcards.POP}_chr{wildcards.CHR}"),
    shell("mv {wildcards.POP}_chr{wildcards.CHR}.anc {output[0]}"),
    shell("mv {wildcards.POP}_chr{wildcards.CHR}.mut {output[1]}")

rule convert_to_tskit:
  input:
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.anc", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}"),
    expand("{tmp_dir}/{POP}/{POP}_chr{CHR}/{POP}_chr{CHR}.mut", tmp_dir=config['tmp_dir'], POP="{POP}", CHR="{CHR}")
  resources:
    time="1:00:00",
    threads=1,
    mem="15gb"
  output:
    expand("{trees_dir}/{POP}/{POP}_chr{CHR}.trees", trees_dir=config['trees_dir'], POP="{POP}", CHR="{CHR}")
  shell:
    "{config[path_to_relate]}/bin/RelateFileFormats --mode ConvertToTreeSequence -i {config[tmp_dir]}/{wildcards.POP}/{wildcards.POP}_chr{wildcards.CHR}/{wildcards.POP}_chr{wildcards.CHR} -o {config[trees_dir]}/{wildcards.POP}/{wildcards.POP}_chr{wildcards.CHR}"

