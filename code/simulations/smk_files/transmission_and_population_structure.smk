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
  [0.75, 0.0, 0.25],
  [0.5, 0.0, 0.5],
  [0.0, 0.0, 1.0]
  ]

tm_K = [1000, 5000, 10000]

tm_len = [1e7]

tm_SS = [100]

tm_seeds = [str(i) for i in range(0,32)]

rule tm_df:
  input:
    expand("{results_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}.results", results_dir=config['results_dir'], tr = tm_tr, K = tm_K, SS = tm_SS, LEN = tm_len)
  resources:
    time="0:30:00",
    threads=1,
    mem="3gb"
  output:
    expand("{results_dir}/dataframes/transmission.pkl", results_dir=config['results_dir'])
  shell:
    '''
    python scripts/get_df.py -output_file {output} -input_files {input}
    '''

rule calculate_tree_distances_neutral_arg:
  input:
    expand("{sim_dir}/neutral_arg/host_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/neutral_arg/symb_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/neutral_arg/pairs_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  params:
  resources:
    time="48:00:00",
    threads=1,
    mem="5gb"
  output:
    expand("{results_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.results", results_dir=config['results_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  shell:
    '''
    python scripts/arg_tree_comparison.py --host_tree {input[0]} --symb_tree {input[1]} --pairs {input[2]} --output_file {output} --sample_size {wildcards.SS} --stratify_col "subpop"
    '''

rule aggregate_tree_distances_neutral_arg:
  input:
    expand("{results_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{seeds}.results", results_dir=config['results_dir'], seeds=tm_seeds, PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  params:
  resources:
    time="00:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}.results", results_dir=config['results_dir'], ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  shell:
    '''
    cat {input} > {output}
    '''

rule simulate_neutral_arg:
  input:
  resources:
    time="12:00:00",
    threads=1,
    mem="1gb",
  output:
    expand("{sim_dir}/neutral_arg/host_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/neutral_arg/symb_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/neutral_arg/pairs_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  shell:
    "slim -s {wildcards.SEED} "
    "-d \"host_tree_outfile='{output[0]}'\" "
    "-d \"symb_tree_outfile='{output[1]}'\" "
    "-d \"pairs_outfile='{output[2]}'\" "
    "-d K={wildcards.K} "
    "-d LEN={wildcards.LEN} "
    "-d p_mother={wildcards.PM} "
    "-d p_father={wildcards.PF} "
    "-d p_random={wildcards.PR} "
    "scripts/slim/neutral_arg_nonWF.slim"

##########################################
### Transmission mode with filter
##########################################

tm_w_filter_tr = [
  [1.0, 0.0, 0.0],
  ]

tm_w_filter_K = [1000, 5000, 10000]

tm_w_filter_len = [1e7]

tm_w_filter_filter = ['parent', 'grandparent']

tm_w_filter_SS = [100]

tm_w_filter_seeds = [str(i) for i in range(0,32)]

rule tm_w_filter_df:
  input:
    expand("{results_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}_FILTER_{FILTER}.results", results_dir=config['results_dir'], tr = tm_w_filter_tr, SS = tm_w_filter_SS, K = tm_w_filter_K, LEN = tm_w_filter_len, FILTER = tm_w_filter_filter)
  resources:
    time="0:30:00",
    threads=1,
    mem="3gb"
  output:
    expand("{results_dir}/dataframes/transmission_w_filter.pkl", results_dir=config['results_dir'])
  shell:
    '''
    python scripts/get_df.py -output_file {output} -input_files {input}
    '''

rule calculate_tree_distances_neutral_arg_w_filter:
  input:
    expand("{sim_dir}/neutral_arg/host_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/neutral_arg/symb_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/neutral_arg/pairs_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  params:
  resources:
    time="48:00:00",
    threads=1,
    mem="5gb"
  output:
    expand("{results_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_FILTER_{FILTER}_seed_{SEED}.results", results_dir=config['results_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}", FILTER="{FILTER}")
  shell:
    '''
    python scripts/arg_tree_comparison.py --host_tree {input[0]} --symb_tree {input[1]} --pairs {input[2]} --filter {wildcards.FILTER} --output_file {output} --sample_size {wildcards.SS} --stratify_col "subpop"
    '''

rule aggregate_tree_distances_neutral_arg_w_filter:
  input:
    expand("{results_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_FILTER_{FILTER}_seed_{SEED}.results", results_dir=config['results_dir'], SEED=tm_w_filter_seeds, PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}", FILTER="{FILTER}")
  params:
  resources:
    time="00:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/neutral_arg/K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_FILTER_{FILTER}.results", results_dir=config['results_dir'], ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}", FILTER="{FILTER}")
  shell:
    '''
    cat {input} > {output}
    '''

##########################################
### Migration
##########################################

mig_tr = [
  [1.0, 0.0, 0.0],
  [0.0, 0.0, 1.0]
  ]

mig_K = [5000]

mig_M = ["0.0005", "0.001", "0.005", "0.01"]

mig_len = [1e7]

mig_SS = [200]

mig_seeds = [str(i) for i in range(0,32)]

rule mig_df:
  input:
    expand("{results_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}.results", results_dir=config['results_dir'], tr = mig_tr, K = mig_K, SS = mig_SS, M = mig_M, LEN = mig_len)
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/dataframes/migration.pkl", results_dir=config['results_dir'])
  shell:
    '''
    python scripts/get_df.py -output_file {output} -input_files {input}
    '''

rule calculate_tree_distances_migration_arg:
  input:
    expand("{sim_dir}/migration_arg/host_tree_M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/migration_arg/symb_tree_M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/migration_arg/pairs_M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  params:
  resources:
    time="48:00:00",
    threads=1,
    mem="5gb"
  output:
    expand("{results_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.results", results_dir=config['results_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  shell:
    '''
    python scripts/arg_tree_comparison.py --host_tree {input[0]} --symb_tree {input[1]} --pairs {input[2]} --output_file {output} --sample_size {wildcards.SS} --stratify_col "subpop"
    '''

rule aggregate_tree_distances_migration_arg:
  input:
    expand("{results_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{seeds}.results", results_dir=config['results_dir'], seeds=mig_seeds, M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  params:
  resources:
    time="00:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}.results", results_dir=config['results_dir'], M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  shell:
    '''
    cat {input} > {output}
    '''

rule simulate_migration_arg:
  input:
  resources:
    time="24:00:00",
    threads=1,
    mem="1gb"
  output:
    expand("{sim_dir}/migration_arg/host_tree_M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/migration_arg/symb_tree_M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/migration_arg/pairs_M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  shell:
    "slim -s {wildcards.SEED} "
    "-d \"host_tree_outfile='{output[0]}'\" "
    "-d \"symb_tree_outfile='{output[1]}'\" "
    "-d \"pairs_outfile='{output[2]}'\" "
    "-d K={wildcards.K} "
    "-d LEN={wildcards.LEN} "
    "-d p_mother={wildcards.PM} "
    "-d p_father={wildcards.PF} "
    "-d p_random={wildcards.PR} "
    "-d migration_rate={wildcards.M} "
    "scripts/slim/migration_arg_nonWF.slim"


##########################################
### Sympatric Migration
##########################################

symp_mig_tr = [
  [0.0, 0.0, 1.0]
  ]

symp_mig_K = [5000]

symp_mig_M = ["0.0005", "0.001", "0.005", "0.01"]

symp_mig_len = [1e7]

symp_mig_SS = [200]

symp_mig_seeds = [str(i) for i in range(0,32)]

rule symp_mig_df:
  input:
    expand("{results_dir}/symp_migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}.results", results_dir=config['results_dir'], tr = symp_mig_tr, K = symp_mig_K, SS = symp_mig_SS, M = symp_mig_M, LEN = symp_mig_len)
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/dataframes/symp_migration.pkl", results_dir=config['results_dir'])
  shell:
    '''
    python scripts/get_df.py -output_file {output} -input_files {input}
    '''

rule calculate_tree_distances_symp_migration_arg:
  input:
    expand("{sim_dir}/symp_migration_arg/host_tree_M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/symp_migration_arg/symb_tree_M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/symp_migration_arg/pairs_M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  params:
  resources:
    time="48:00:00",
    threads=1,
    mem="5gb"
  output:
    expand("{results_dir}/symp_migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.results", results_dir=config['results_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  shell:
    '''
    python scripts/arg_tree_comparison.py --host_tree {input[0]} --symb_tree {input[1]} --pairs {input[2]} --output_file {output} --sample_size {wildcards.SS} --stratify_col "subpop"
    '''

rule aggregate_tree_distances_symp_migration_arg:
  input:
    expand("{results_dir}/symp_migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{seeds}.results", results_dir=config['results_dir'], seeds=symp_mig_seeds, M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  params:
  resources:
    time="00:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/symp_migration_arg/M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}.results", results_dir=config['results_dir'], M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  shell:
    '''
    cat {input} > {output}
    '''

rule simulate_symp_migration_arg:
  input:
  resources:
    time="24:00:00",
    threads=1,
    mem="1gb"
  output:
    expand("{sim_dir}/symp_migration_arg/host_tree_M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/symp_migration_arg/symb_tree_M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/symp_migration_arg/pairs_M_{M}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", M="{M}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  shell:
    "slim -s {wildcards.SEED} "
    "-d \"host_tree_outfile='{output[0]}'\" "
    "-d \"symb_tree_outfile='{output[1]}'\" "
    "-d \"pairs_outfile='{output[2]}'\" "
    "-d K={wildcards.K} "
    "-d LEN={wildcards.LEN} "
    "-d p_mother={wildcards.PM} "
    "-d p_father={wildcards.PF} "
    "-d p_random={wildcards.PR} "
    "-d migration_rate={wildcards.M} "
    "scripts/slim/sympatric_migration_arg_nonWF.slim"

##########################################
### Pop Split ARG
##########################################

ps_tr = [
  [1.0, 0.0, 0.0],
  [0.0, 0.0, 1.0]
  ]

ps_K = [5000]

ps_len = [1e7]

ps_sample_times = ["1", "10", "100", "1000"]

ps_SS = [200]

ps_seeds = [str(i) for i in range(0,32)]

rule pop_split_arg:
  input:
    expand("{results_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}.results", results_dir=config['results_dir'], tr = ps_tr, K = ps_K, SS = ps_SS, LEN = ps_len, ST = ps_sample_times)
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/dataframes/population_split.pkl", results_dir=config['results_dir'])
  shell:
    '''
    python scripts/get_df.py -output_file {output} -input_files {input}
    '''

rule calculate_tree_distances_population_split_arg:
  input:
    expand("{sim_dir}/population_split_arg/host_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/population_split_arg/symb_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/population_split_arg/pairs_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  params:
  resources:
    time="48:00:00",
    threads=1,
    mem="5gb"
  output:
    expand("{results_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.results", results_dir=config['results_dir'], SEED ="{SEED}", ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  shell:
    '''
    python scripts/arg_tree_comparison.py --host_tree {input[0]} --symb_tree {input[1]} --pairs {input[2]} --output_file {output} --sample_size {wildcards.SS} --stratify_col "subpop"
    '''

rule aggregate_tree_distances_population_split_arg:
  input:
    expand("{results_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{seeds}.results", results_dir=config['results_dir'], seeds=ps_seeds, ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  params:
  resources:
    time="00:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}.results", results_dir=config['results_dir'], ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  shell:
    '''
    cat {input} > {output}
    '''

rule simulate_population_split_arg:
  input:
  resources:
    time="24:00:00",
    threads=1,
    mem="1gb"
  params:
    expand("{sim_dir}/population_split_arg/host_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/population_split_arg/symb_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/population_split_arg/pairs_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  output:
    expand("{sim_dir}/population_split_arg/host_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}", ST=ps_sample_times, PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/population_split_arg/symb_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}", ST=ps_sample_times, PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/population_split_arg/pairs_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", ST=ps_sample_times, PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  shell:
    "slim -s {wildcards.SEED} "
    "-d \"host_tree_header='{params[0]}'\" "
    "-d \"symb_tree_header='{params[1]}'\" "
    "-d \"pairs_header='{params[2]}'\" "
    "-d K={wildcards.K} "
    "-d LEN={wildcards.LEN} "
    "-d p_mother={wildcards.PM} "
    "-d p_father={wildcards.PF} "
    "-d p_random={wildcards.PR} "
    "-d sample_time=1000.0 "
    "scripts/slim/population_split_arg_nonWF.slim"

##########################################
### Sympatric Pop Split ARG
##########################################

symp_ps_tr = [
  [0.0, 0.0, 1.0]
  ]

symp_ps_K = [5000]

symp_ps_len = [1e7]

symp_ps_sample_times = ["1", "10", "100", "1000"]

symp_ps_SS = [200]

symp_ps_seeds = [str(i) for i in range(0,32)]

rule symp_pop_split_arg:
  input:
    expand("{results_dir}/symp_population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}.results", results_dir=config['results_dir'], tr = symp_ps_tr, K = symp_ps_K, SS = symp_ps_SS, LEN = symp_ps_len, ST = symp_ps_sample_times)
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/dataframes/symp_population_split.pkl", results_dir=config['results_dir'])
  shell:
    '''
    python scripts/get_df.py -output_file {output} -input_files {input}
    '''

rule calculate_tree_distances_symp_population_split_arg:
  input:
    expand("{sim_dir}/symp_population_split_arg/host_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/symp_population_split_arg/symb_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/symp_population_split_arg/pairs_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  params:
  resources:
    time="48:00:00",
    threads=1,
    mem="5gb"
  output:
    expand("{results_dir}/symp_population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.results", results_dir=config['results_dir'], SEED ="{SEED}", ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  shell:
    '''
    python scripts/arg_tree_comparison.py --host_tree {input[0]} --symb_tree {input[1]} --pairs {input[2]} --output_file {output} --sample_size {wildcards.SS} --stratify_col "subpop"
    '''

rule aggregate_tree_distances_symp_population_split_arg:
  input:
    expand("{results_dir}/symp_population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{seeds}.results", results_dir=config['results_dir'], seeds=symp_ps_seeds, ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  params:
  resources:
    time="00:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/symp_population_split_arg/ST_{ST}_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}.results", results_dir=config['results_dir'], ST="{ST}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}")
  shell:
    '''
    cat {input} > {output}
    '''

rule simulate_symp_population_split_arg:
  input:
  resources:
    time="24:00:00",
    threads=1,
    mem="1gb"
  params:
    expand("{sim_dir}/symp_population_split_arg/host_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/symp_population_split_arg/symb_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/symp_population_split_arg/pairs_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_", sim_dir=config['sim_dir'], SEED ="{SEED}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  output:
    expand("{sim_dir}/symp_population_split_arg/host_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}", ST=symp_ps_sample_times, PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/symp_population_split_arg/symb_tree_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}", ST=symp_ps_sample_times, PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
    expand("{sim_dir}/symp_population_split_arg/pairs_K_{K}_SS_{SS}_LEN_{LEN}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", ST=symp_ps_sample_times, PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}", LEN="{LEN}"),
  shell:
    "slim -s {wildcards.SEED} "
    "-d \"host_tree_header='{params[0]}'\" "
    "-d \"symb_tree_header='{params[1]}'\" "
    "-d \"pairs_header='{params[2]}'\" "
    "-d K={wildcards.K} "
    "-d LEN={wildcards.LEN} "
    "-d p_mother={wildcards.PM} "
    "-d p_father={wildcards.PF} "
    "-d p_random={wildcards.PR} "
    "-d sample_time=1000.0 "
    "scripts/slim/sympatric_population_split_arg_nonWF.slim"