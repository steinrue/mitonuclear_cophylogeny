configfile: "config.yaml"

wildcard_constraints:
  TIPS = "|".join(["50", "100", "200"]),
  PERMS = "|".join(["100000"]),
  IDX = "|".join([str(i) for i in range(0,1000)])

num_tips = ["50", "100", "200"]
num_perms = ["100000"]

idx_to_run = [str(i) for i in range(0,64)]

rule sim_binary_trees:
    input:
    resources:
        time="1:00:00",
        threads=1,
        mem="1gb"
    output:
        expand("{sim_dir}/simulated_binary_trees_w_{TIPS}_tips.list", sim_dir=config['sim_dir'], TIPS="{TIPS}")
    shell:
        '''
        python3 scripts/sim_binary_trees.py -o {output} --num_trees 1000 --num_leaves {wildcards.TIPS}
        '''


rule aggregate_ind_tree_results:
  input:
    expand("{tmp_dir}/cophylo_pvals/tips_{TIPS}_perms_{PERMS}_focal_{IDX}.results", tmp_dir=config['tmp_dir'], TIPS="{TIPS}", PERMS ="{PERMS}", IDX = idx_to_run)
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{tmp_dir}/cophylo_pvals/tips_{TIPS}_perms_{PERMS}.results", tmp_dir=config['tmp_dir'], TIPS="{TIPS}", PERMS ="{PERMS}")
  shell:
    '''
    cat {input} > {output}
    '''

rule get_tree_distances:
  input:
    expand("{sim_dir}/simulated_binary_trees_w_{TIPS}_tips.list", sim_dir=config['sim_dir'], TIPS="{TIPS}")
  resources:
    time="5:30:00",
    threads=1,
    mem="3gb"
  output:
    expand("{tmp_dir}/cophylo_pvals/tips_{TIPS}_perms_{PERMS}_focal_{IDX}.results", tmp_dir=config['tmp_dir'], TIPS="{TIPS}", PERMS ="{PERMS}", IDX = "{IDX}"),
  shell:
    '''
    python3 scripts/tree_comp.py -tree_list {input[0]} -focal_tree_idx {wildcards.IDX} -num_perms {wildcards.PERMS} -o {output[0]}    
    '''

rule fit_nig_pvalues:
  input:
    expand("{tmp_dir}/cophylo_pvals/tips_{TIPS}_perms_{PERMS}.results", tmp_dir=config['tmp_dir'], TIPS="{TIPS}", PERMS ="{PERMS}")
  resources:
    time="4:30:00",
    threads=1,
    mem="5gb"
  output:
    expand("{results_dir}/tips_{TIPS}_perms_{PERMS}.pkl", results_dir=config['results_dir'], TIPS="{TIPS}", PERMS ="{PERMS}")
  shell:
    '''
    python3 scripts/norminvgauss.py --input {input} --output {output}
    '''