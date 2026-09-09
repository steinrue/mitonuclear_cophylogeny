configfile: "config.yaml"

wildcard_constraints:
    PM = r"\d+(\.\d+)?",
    PF = r"\d+(\.\d+)?",
    PR = r"\d+(\.\d+)?",
    FILTER = "parent|grandparent",
    SEED = r"\d+"

##########################################
### Transmission mode
##########################################

tm_tr = [
  [1.0, 0.0, 0.0],
  ]

tm_K = [1000, 5000, 10000]

tm_len = [1e7]

tm_SS = [100]

tm_seeds = [str(i) for i in range(0, 32)]

rule tm_fis_df:
  input:
    expand("{sim_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}.fis", sim_dir=config['sim_dir'], tr = tm_tr, K = tm_K, SS = tm_SS, LEN = tm_len),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/dataframes/fis_tm.pkl", results_dir=config['results_dir'])
  shell:
    '''
    python scripts/get_df.py -output_file {output} -input_files {input}
    '''

rule get_vcf_neutral:
  input:
    expand("{sim_dir}/neutral_arg/host_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/neutral_arg/pairs_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  params:
  resources:
    time="12:00:00",
    threads=1,
    mem="3gb"
  output:
    expand("{sim_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.vcf", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  shell:
    '''
    python scripts/get_vcf.py --host_tree {input[0]} --pairs {input[1]} --output_file {output} --sample_size {wildcards.SS} --stratify_col "subpop"
    '''

rule get_biallelic_vcf_neutral:
  input:
    expand("{sim_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.vcf", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  resources:
    time="2:00:00",
    threads=1,
    mem="10gb"
  output:
    expand("{sim_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_filtered.vcf", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  shell:
    "bcftools view -m2 -M2 -v snps --min-af 0.000001 {input} -Oz -o {output} "

rule get_hardy_neutral:
  input:
    expand("{sim_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_filtered.vcf", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  params:
    outfile = expand("{sim_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  resources:
    time="1:00:00",
    threads=1,
    mem="12gb"
  output:
    expand("{sim_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.hardy", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  shell:
    "plink2 --vcf {input} --vcf-half-call h "
    "--hardy midp "
    "--out {params.outfile} "

rule get_fis_neutral:
  input:
    expand("{sim_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.hardy", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  resources:
    time="1:00:00",
    threads=1,
    mem="24gb"
  output:
    expand("{sim_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.fis", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  shell:
    "python scripts/fis_jackknife.py {input} > {output}"

rule get_final_fis_neutral:
  input:
    expand("{sim_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.fis", sim_dir=config['sim_dir'], SEED = tm_seeds, PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  resources:
    time="5:00:00",
    threads=1,
    mem="24gb"
  output:
    expand("{sim_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}.fis", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  shell:
    """
    (echo -e "population\tfis\tlower\tupper"; cat {input}) > {output}
    """

##########################################
### Migration
##########################################

mig_tr = [
  [1.0, 0.0, 0.0],
  ]

mig_K = [5000]

mig_M = ["0.0005", "0.001", "0.005", "0.01"]

mig_len = [1e7]

mig_SS = [200]

mig_seeds = [str(i) for i in range(0, 32)]

rule mig_fis_df:
  input:
    expand("{sim_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}.fis", sim_dir=config['sim_dir'], tr = mig_tr, K = mig_K, SS = mig_SS, M = mig_M, LEN = mig_len),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/dataframes/fis_mig.pkl", results_dir=config['results_dir'])
  shell:
    '''
    python scripts/get_df.py -output_file {output} -input_files {input}
    '''

rule get_vcf_migration:
  input:
    expand("{sim_dir}/migration_arg/host_tree_M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/migration_arg/pairs_M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  params:
  resources:
    time="12:00:00",
    threads=1,
    mem="3gb"
  output:
    expand("{sim_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.vcf", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  shell:
    '''
    python scripts/get_vcf.py --host_tree {input[0]} --pairs {input[1]} --output_file {output} --sample_size {wildcards.SS} --stratify_col "subpop"
    '''

rule get_biallelic_vcf_migration:
  input:
    expand("{sim_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.vcf", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  resources:
    time="2:00:00",
    threads=1,
    mem="10gb"
  output:
    expand("{sim_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_filtered.vcf.gz", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  shell:
    "bcftools view -m2 -M2 -v snps --min-af 0.000001 {input} -Oz -o {output} "


rule get_hardy_migration:
  input:
    expand("{sim_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_filtered.vcf.gz", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  params:
    outfile = expand("{sim_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  resources:
    time="1:00:00",
    threads=1,
    mem="12gb"
  output:
    expand("{sim_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.hardy", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  shell:
    "plink2 --vcf {input} --vcf-half-call h "
    "--hardy midp "
    "--out {params.outfile} "

rule get_fis_migration:
  input:
    expand("{sim_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.hardy", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  resources:
    time="1:00:00",
    threads=1,
    mem="24gb"
  output:
    expand("{sim_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.fis", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  shell:
    "python scripts/fis_jackknife.py {input} > {output}"

rule get_final_fis_migration:
  input:
    expand("{sim_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.fis", sim_dir=config['sim_dir'], SEED =mig_seeds, M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  resources:
    time="5:00:00",
    threads=1,
    mem="24gb"
  output:
    expand("{sim_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}.fis", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  shell:
    """
    (echo -e "population\tfis\tlower\tupper"; cat {input}) > {output}
    """

##########################################
### Pop Split ARG
##########################################

ps_tr = [
  [1.0, 0.0, 0.0]
  ]

ps_K = [5000]

ps_len = [1e6]

ps_sample_times = ["1", "10", "100", "1000"]

ps_SS = [200]

ps_seeds = [str(i) for i in range(0, 32)]

rule population_split_fis_df:
  input:
    expand("{sim_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}.fis", sim_dir=config['sim_dir'], tr = ps_tr, K = ps_K, SS = ps_SS, LEN = ps_len, ST = ps_sample_times)
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/dataframes/fis_pop_split.pkl", results_dir=config['results_dir'])
  shell:
    '''
    python scripts/get_df.py -output_file {output} -input_files {input}
    '''

rule get_vcf_population_split:
  input:
    expand("{sim_dir}/population_split_arg/host_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}", ST = "{ST}"),
    expand("{sim_dir}/population_split_arg/pairs_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}", ST = "{ST}")
  params:
  resources:
    time="12:00:00",
    threads=1,
    mem="3gb"
  output:
    expand("{sim_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.vcf", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}", ST = "{ST}")
  shell:
    '''
    python scripts/get_vcf.py --host_tree {input[0]} --pairs {input[1]} --output_file {output} --sample_size {wildcards.SS} --stratify_col "subpop"
    '''

rule get_biallelic_vcf_population_split:
  input:
    expand("{sim_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.vcf", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}", ST = "{ST}")
  resources:
    time="2:00:00",
    threads=1,
    mem="10gb"
  output:
    expand("{sim_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_filtered.vcf.gz", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}", ST = "{ST}")
  shell:
    "bcftools view -m2 -M2 -v snps --min-af 0.000001 {input} -Oz -o {output} "

rule get_hardy_population_split:
  input:
    expand("{sim_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_filtered.vcf.gz", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}", ST = "{ST}")
  params:
    outfile = expand("{sim_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}", ST = "{ST}"),
  resources:
    time="1:00:00",
    threads=1,
    mem="12gb"
  output:
    expand("{sim_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.hardy", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}", ST = "{ST}"),
  shell:
    "plink2 --vcf {input} --vcf-half-call h "
    "--hardy midp "
    "--out {params.outfile} "

rule get_fis_population_split:
  input:
    expand("{sim_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.hardy", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}", ST = "{ST}"),
  resources:
    time="1:00:00",
    threads=1,
    mem="24gb"
  output:
    expand("{sim_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.fis", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  shell:
    "python scripts/fis_jackknife.py {input} > {output}"

rule get_final_fis_population_split:
  input:
    expand("{sim_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.fis", sim_dir=config['sim_dir'], SEED = ps_seeds, ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  resources:
    time="5:00:00",
    threads=1,
    mem="24gb"
  output:
    expand("{sim_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}.fis", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  shell:
    """
    (echo -e "population\tfis\tlower\tupper"; cat {input}) > {output}
    """

##########################################
### Admixture
##########################################

f_adm_tr = [
  [1.0, 0.0, 0.0]
  ]

f_adm_S = [
  [0.00, 0.00],
  ]

f_adm_K = [5000]

f_adm_len = [1e7]

f_adm_mut = [1e-8]

f_adm_sample_times = ["1", "2", "3", "5", "10", "25", "100"]

f_adm_SS = [100]

f_adm_prop = [.5]

f_adm_seeds = [str(i) for i in range(0, 32)]

rule admixture_fis_df:
  input:
    expand("{sim_dir}/admixture/S1_{S[0]}_S2_{S[1]}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}_ST_{ST}.fis", sim_dir=config['sim_dir'],MUT = f_adm_mut, S = f_adm_S, ST=f_adm_sample_times, tr = f_adm_tr, K = f_adm_K, SS = f_adm_SS, LEN = f_adm_len, PROP = f_adm_prop),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/dataframes/fis_admixture.pkl", results_dir=config['results_dir'])
  shell:
    '''
    python scripts/get_df.py -output_file {output} -input_files {input}
    '''

rule get_vcf_admixture:
  input:
    expand("{sim_dir}/admixture/host_tree_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}"),
    expand("{sim_dir}/admixture/pairs_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}"),
  params:
  resources:
    time="12:00:00",
    threads=1,
    mem="3gb"
  output:
    expand("{sim_dir}/admixture/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.vcf", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}", SS = "{SS}"),
  shell:
    '''
    python scripts/get_vcf.py --host_tree {input[0]} --pairs {input[1]} --output_file {output} --sample_size {wildcards.SS} --stratify_col "subpop"
    '''

rule get_biallelic_vcf_admixture:
  input:
    expand("{sim_dir}/admixture/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.vcf", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}", SS = "{SS}"),
  resources:
    time="2:00:00",
    threads=1,
    mem="10gb"
  output:
    expand("{sim_dir}/admixture/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}_filtered.vcf.gz", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}", SS = "{SS}"),
  shell:
    "bcftools view -m2 -M2 -v snps --min-af 0.000001 {input} -Oz -o {output} "

rule get_hardy_admixture:
  input:
    expand("{sim_dir}/admixture/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}_filtered.vcf.gz", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}", SS = "{SS}"),
  params:
    outfile = expand("{sim_dir}/admixture/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}", SS = "{SS}"),
  resources:
    time="1:00:00",
    threads=1,
    mem="12gb"
  output:
    expand("{sim_dir}/admixture/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.hardy", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}", SS = "{SS}"),
  shell:
    "plink2 --vcf {input} --vcf-half-call h "
    "--hardy midp "
    "--out {params.outfile} "

rule get_fis_admixture:
  input:
    expand("{sim_dir}/admixture/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.hardy", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}", SS = "{SS}"),
  resources:
    time="1:00:00",
    threads=1,
    mem="24gb"
  output:
    expand("{sim_dir}/admixture/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.fis", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}", SS = "{SS}"),
  shell:
    "python scripts/fis_jackknife.py {input} > {output}"

rule get_final_fis_admixture:
  input:
    expand("{sim_dir}/admixture/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.fis", sim_dir=config['sim_dir'], SEED = f_adm_seeds, ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}", SS = "{SS}"),
  resources:
    time="5:00:00",
    threads=1,
    mem="24gb"
  output:
    expand("{sim_dir}/admixture/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_ST_{ST}.fis", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}", SS = "{SS}"),
  shell:
    """
    (echo -e "population\tfis\tlower\tupper"; cat {input}) > {output}
    """