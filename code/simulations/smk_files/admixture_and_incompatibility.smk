configfile: "config.yaml"

wildcard_constraints:
    PM = r"\d+(\.\d+)?",
    PF = r"\d+(\.\d+)?",
    PR = r"\d+(\.\d+)?",
    SEED = r"\d+"

##########################################
### Allele Incompatability
##########################################

f_ai_tr = [
  [1.0, 0.0, 0.0]
  ]

f_ai_S = [
  [0.0, 1]
  ]

f_ai_K = [10000]

f_ai_len = [1e7]

f_ai_mut = [1e-8]

f_ai_sample_times = ["0", "5", "10", "100"]

f_ai_SS = [100]

f_ai_seeds = [str(i) for i in range(0, 32)]

rule allele_incompatability:
  input:
    expand("{results_dir}/allele_incompatability/S1_{S[0]}_S2_{S[1]}_LEN_{LEN}_MUT_{MUT}_K_{K}_SS_{SS}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}_ST_{ST}.results", results_dir=config['results_dir'], tr = f_ai_tr, K = f_ai_K, ST = f_ai_sample_times, S = f_ai_S, SS = f_ai_SS, LEN = f_ai_len, MUT = f_ai_mut),
    expand("{sim_dir}/allele_incompatability/allele_freqs_S1_{S[0]}_S2_{S[1]}_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED = f_ai_seeds, tr = f_ai_tr, K = f_ai_K, ST = f_ai_sample_times, S = f_ai_S, SS = f_ai_SS, LEN = f_ai_len, MUT = f_ai_mut),
  resources:
    time="0:30:00",
    threads=1,
    mem="500mb"
  output:
    expand("{results_dir}/dataframes/allele_incompatability.pkl", results_dir=config['results_dir'])
  shell:
    '''
    python scripts/get_df.py -output_file {output} -input_files {input} --allele_freq_data
    '''

rule allele_incompatability_aft:
  input:
      expand("{sim_dir}/allele_incompatability/allele_freqs_S1_{S[0]}_S2_{S[1]}_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED = f_ai_seeds, tr = f_ai_tr, K = f_ai_K, ST = f_ai_sample_times, S = f_ai_S, SS = f_ai_SS, LEN = f_ai_len, MUT = f_ai_mut),  
  resources:
    time="0:30:00",
    threads=1,
    mem="3gb"
  output:
    expand("{results_dir}/dataframes/allele_incompatability_afts.pkl", results_dir=config['results_dir'])
  shell:
    '''
    python scripts/get_allele_freqs_df.py --allele-freqs {input} --output-pkl {output}
    '''

rule calculate_tree_distances_allele_incompatability_arg:
  input:
    expand("{sim_dir}/allele_incompatability/host_tree_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
    expand("{sim_dir}/allele_incompatability/symb_tree_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
    expand("{sim_dir}/allele_incompatability/pairs_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
  params:
  resources:
    time="72:00:00",
    threads=1,
    mem="5gb"
  output:
    expand("{results_dir}/allele_incompatability/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.results", results_dir=config['results_dir'], SEED ="{SEED}", ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", ITR = "{ITR}", SS="{SS}")
  shell:
    '''
    python scripts/arg_tree_comparison.py --host_tree {input[0]} --symb_tree {input[1]} --pairs {input[2]} --output_file {output} --sample_size {wildcards.SS} --stratify_col "allele"
    '''

rule aggregate_tree_distances_allele_incompatability:
  input:
    expand("{results_dir}/allele_incompatability/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{seeds}_ST_{ST}.results", results_dir=config['results_dir'], seeds=f_ai_seeds, ST="{ST}", S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}")
  params:
  resources:
    time="0:10:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/allele_incompatability/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_ST_{ST}.results", results_dir=config['results_dir'], ST="{ST}", S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}")
  shell:
    '''
    cat {input} > {output}
    '''

rule simulate_allele_incompatability_arg:
  input:
    expand("{sim_dir}/selection/host_savepoint_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
    expand("{sim_dir}/selection/symb_savepoint_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
  params:
    expand("{sim_dir}/allele_incompatability/host_tree_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST=f_ai_sample_times,  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
    expand("{sim_dir}/allele_incompatability/symb_tree_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST=f_ai_sample_times,  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
    expand("{sim_dir}/allele_incompatability/pairs_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_", sim_dir=config['sim_dir'], SEED ="{SEED}", ST=f_ai_sample_times,  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
  resources:
    time="8:00:00",
    threads=1,
    mem="12gb"
  output:
    expand("{sim_dir}/allele_incompatability/allele_freqs_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
    expand("{sim_dir}/allele_incompatability/host_tree_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST=f_ai_sample_times,  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
    expand("{sim_dir}/allele_incompatability/symb_tree_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST=f_ai_sample_times,  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
    expand("{sim_dir}/allele_incompatability/pairs_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", ST=f_ai_sample_times,  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),

  shell:
    "slim -s {wildcards.SEED} "
    "-d \"host_savepoint='{input[0]}'\" "
    "-d \"symb_savepoint='{input[1]}'\" "
    "-d \"host_tree_header='{params[0]}'\" "
    "-d \"symb_tree_header='{params[1]}'\" "
    "-d \"pairs_header='{params[2]}'\" "
    "-d \"allele_freq_outfile='{output[0]}'\" "
    "-d K={wildcards.K} "
    "-d p_mother={wildcards.PM} "
    "-d p_father={wildcards.PF} "
    "-d p_random={wildcards.PR} "
    "-d MUT={wildcards.MUT} "
    "-d LEN={wildcards.LEN} "
    "-d S1={wildcards.S1} "
    "-d S2={wildcards.S2} "
    "scripts/slim/allele_incompatability_arg_nonWF.slim"

rule simulate_savepoint:
  input:
  resources:
    time="48:00:00",
    threads=1,
    mem="12gb"
  output:
    expand("{sim_dir}/selection/host_savepoint_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
    expand("{sim_dir}/selection/symb_savepoint_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
    expand("{sim_dir}/selection/host_mutations_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}",  MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
    expand("{sim_dir}/selection/symb_mutations_LEN_{LEN}_MUT_{MUT}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}",  MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}"),
  shell:
    "slim -s {wildcards.SEED} "
    "-d \"host_savepoint='{output[0]}'\" "
    "-d \"symb_savepoint='{output[1]}'\" "
    "-d \"host_mut_freqs='{output[2]}'\" "
    "-d \"symb_mut_freqs='{output[3]}'\" "
    "-d K={wildcards.K} "
    "-d p_mother={wildcards.PM} "
    "-d p_father={wildcards.PF} "
    "-d p_random={wildcards.PR} "
    "-d MUT={wildcards.MUT} "
    "-d LEN={wildcards.LEN} "
    "scripts/slim/neutral_mutations_arg_nonWF.slim"

##########################################
### Admixture
##########################################

f_adm_tr = [
  [1.0, 0.0, 0.0]
  ]

f_adm_S = [
  [0.00, 0.00],
  [0.0, 1],
  ]

f_adm_K = [5000]

f_adm_len = [1e7]

f_adm_mut = [0.0]

f_adm_sample_times = ["1", "2", "3", "5", "10", "25", "100"]

f_adm_SS = [100]

f_adm_prop = [.1, .2, .5]

f_adm_seeds = [str(i) for i in range(0, 32)]

rule admixture:
  input:
    expand("{results_dir}/admixture/S1_{S[0]}_S2_{S[1]}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}_ST_{ST}.results", results_dir=config['results_dir'], tr = f_adm_tr, K = f_adm_K, ST = f_adm_sample_times, S = f_adm_S, SS = f_adm_SS, LEN = f_adm_len, MUT = f_adm_mut, PROP = f_adm_prop),
    expand("{sim_dir}/admixture/allele_freqs_S1_{S[0]}_S2_{S[1]}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED =f_ai_seeds, tr = f_adm_tr, K = f_adm_K, ST = f_adm_sample_times, S = f_adm_S, SS = f_adm_SS, LEN = f_adm_len, MUT = f_adm_mut, PROP = f_adm_prop),
  resources:
    time="0:30:00",
    threads=1,
    mem="3gb"
  output:
    expand("{results_dir}/dataframes/admixture.pkl", results_dir=config['results_dir'])
  shell:
    '''
    python scripts/get_df.py -output_file {output} -input_files {input} --allele_freq_data
    '''

rule admixture_aft:
  input:
      expand("{sim_dir}/admixture/allele_freqs_S1_{S[0]}_S2_{S[1]}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_PM_{tr[0]}_PF_{tr[1]}_PR_{tr[2]}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED =f_ai_seeds, tr = f_adm_tr, K = f_adm_K, ST = f_adm_sample_times, S = f_adm_S, SS = f_adm_SS, LEN = f_adm_len, MUT = f_adm_mut, PROP = f_adm_prop),
  resources:
    time="0:30:00",
    threads=1,
    mem="3gb"
  output:
    expand("{results_dir}/dataframes/admixture_afts.pkl", results_dir=config['results_dir'])
  shell:
    '''
    python scripts/get_allele_freqs_df.py --allele-freqs {input} --output-pkl {output}
    '''

rule calculate_tree_distances_admixture_arg:
  input:
    expand("{sim_dir}/admixture/host_tree_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}"),
    expand("{sim_dir}/admixture/symb_tree_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}"),
    expand("{sim_dir}/admixture/pairs_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}"),
  params:
  resources:
    time="72:00:00",
    threads=1,
    mem="5gb"
  output:
    expand("{results_dir}/admixture/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.results", results_dir=config['results_dir'], SEED ="{SEED}", ST="{ST}",  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}", SS="{SS}")
  shell:
    '''
    python scripts/arg_tree_comparison.py --host_tree {input[0]} --symb_tree {input[1]} --pairs {input[2]} --output_file {output} --sample_size {wildcards.SS} --stratify_col "subpop" --recapitate
    '''

rule aggregate_tree_distances_admixture:
  input:
    expand("{results_dir}/admixture/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{seeds}_ST_{ST}.results", results_dir=config['results_dir'], seeds=f_adm_seeds, ST="{ST}", S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PROP = "{PROP}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", SS="{SS}")
  params:
  resources:
    time="0:10:00",
    threads=1,
    mem="1gb"
  output:
    expand("{results_dir}/admixture/S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_SS_{SS}_PM_{PM}_PF_{PF}_PR_{PR}_ST_{ST}.results", results_dir=config['results_dir'], ST="{ST}", S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}", SS="{SS}")
  shell:
    '''
    cat {input} > {output}
    '''

rule simulate_admixture_arg:
  input:
  params:
    expand("{sim_dir}/admixture/host_tree_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST=f_adm_sample_times,  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}"),
    expand("{sim_dir}/admixture/symb_tree_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST=f_adm_sample_times,  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}"),
    expand("{sim_dir}/admixture/pairs_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_", sim_dir=config['sim_dir'], SEED ="{SEED}", ST=f_adm_sample_times,  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}"),
  resources:
    time="24:00:00",
    threads=1,
    mem="15gb"
  output:
    expand("{sim_dir}/admixture/allele_freqs_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}"),
    expand("{sim_dir}/admixture/host_tree_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST=f_adm_sample_times,  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}"),
    expand("{sim_dir}/admixture/symb_tree_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.trees", sim_dir=config['sim_dir'], SEED ="{SEED}",  ST=f_adm_sample_times,  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}"),
    expand("{sim_dir}/admixture/pairs_S1_{S1}_S2_{S2}_LEN_{LEN}_MUT_{MUT}_PROP_{PROP}_K_{K}_PM_{PM}_PF_{PF}_PR_{PR}_seed_{SEED}_ST_{ST}.txt", sim_dir=config['sim_dir'], SEED ="{SEED}", ST=f_adm_sample_times,  S2="{S2}", S1="{S1}", MUT="{MUT}",  LEN="{LEN}", PM="{PM}", PF="{PF}", PR="{PR}", K="{K}", PROP = "{PROP}"),

  shell:
    "slim -s {wildcards.SEED} "
    "-d \"host_tree_header='{params[0]}'\" "
    "-d \"symb_tree_header='{params[1]}'\" "
    "-d \"pairs_header='{params[2]}'\" "
    "-d \"allele_freq_outfile='{output[0]}'\" "
    "-d K={wildcards.K} "
    "-d p_mother={wildcards.PM} "
    "-d p_father={wildcards.PF} "
    "-d p_random={wildcards.PR} "
    "-d MUT={wildcards.MUT} "
    "-d LEN={wildcards.LEN} "
    "-d S1={wildcards.S1} "
    "-d S2={wildcards.S2} "
    "-d PROP={wildcards.PROP} "
    "scripts/slim/admixture_arg_nonWF.slim"