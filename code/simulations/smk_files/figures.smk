configfile: "config.yaml"

##########################################
### Transmission mode
##########################################

rule tm_bp:
  input:
    expand("{results_dir}/dataframes/transmission_w_filter.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/dataframes/transmission.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/transmission_boxplot.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/transmission_boxplot.py --transmission-w-filter-pkl {input[0]} --transmission-pkl {input[1]} --output-pdf {output}
    '''

rule tm_qq:
  input:
    expand("{results_dir}/dataframes/transmission.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/dataframes/transmission_w_filter.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/transmission_qqplot.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/transmission_qq.py --transmission-pkl {input[0]} --transmission-w-filter-pkl {input[1]} --output-pdf {output}
    '''

rule tm_nf_qq:
  input:
    expand("{results_dir}/dataframes/transmission.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/transmission_no_filter_qqplot.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/transmission_no_filter_qq.py --transmission-pkl {input[0]} --output-pdf {output}
    '''

rule tm_ns_qq:
  input:
    expand("{results_dir}/dataframes/transmission_w_filter.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/transmission_no_siblings_qqplot.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/transmission_no_siblings_qq.py --transmission-w-filter-pkl {input[0]} --output-pdf {output}
    '''

##########################################
### Shared structure
##########################################

rule population_strucutre_bp:
  input:
    expand("{results_dir}/dataframes/migration.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/dataframes/symp_migration.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/dataframes/population_split.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/dataframes/symp_population_split.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/dataframes/admixture.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/population_structure_boxplot.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/population_structure_bp.py --migration-pkl {input[0]} --symp-migration-pkl {input[1]} --population-split-pkl {input[2]} --symp-population-split-pkl {input[3]} --admixture-pkl {input[4]} --output-pdf {output}
    '''

rule population_strucutre_qq:
  input:
    expand("{results_dir}/dataframes/migration.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/dataframes/population_split.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/dataframes/admixture.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/population_structure_qq.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/population_structure_qq.py --migration-pkl {input[0]} --population-split-pkl {input[1]} --admixture-pkl {input[2]} --output-pdf {output}
    '''

rule migration_bp:
  input:
    expand("{results_dir}/dataframes/migration.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/dataframes/symp_migration.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/migration_boxplot.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/migration_boxplot.py --migration-pkl {input[0]} --symp-migration-pkl {input[1]} --output-pdf {output}
    '''

rule migration_qq:
  input:
    expand("{results_dir}/dataframes/migration.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/migration_qqplot.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/migration_qq.py --migration-pkl {input[0]} --output-pdf {output}
    '''

rule population_split_bp:
  input:
    expand("{results_dir}/dataframes/population_split.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/dataframes/symp_population_split.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/population_split_boxplot.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/population_split_boxplot.py --population-split-pkl {input[0]} --symp-population-split-pkl {input[1]} --output-pdf {output}
    '''

rule population_split_qq:
  input:
    expand("{results_dir}/dataframes/population_split.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/population_split_qqplot.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/population_split_qq.py --population-split-pkl {input[0]} --output-pdf {output}
    '''

rule admixture_bp:
  input:
    expand("{results_dir}/dataframes/admixture.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/admixture_boxplot.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/admixture_boxplot.py --admixture-pkl {input[0]} --output-pdf {output}
    '''

rule admixture_qq:
  input:
    expand("{results_dir}/dataframes/admixture.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/admixture_qqplot.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/admixture_qq.py --admixture-pkl {input[0]} --output-pdf {output}
    '''

##########################################
### Allele incompatibility
##########################################

rule allele_incompatability_aft_fig:
  input:
    expand("{results_dir}/dataframes/allele_incompatability_afts.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/allele_incompatability_aft.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/allele_incompatibility_aft.py --allele-incompatability-afts-pkl {input} --output-pdf {output}
    '''

rule allele_incompatability_manhat:
  input:
    expand("{results_dir}/dataframes/allele_incompatability_afts.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/dataframes/allele_incompatability.pkl", results_dir=config['results_dir'])
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/allele_incompatability_manhattan.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/allele_incompatability_manhattan.py --allele-incompatability-afts-pkl {input[0]} --allele-incompatability-pkl {input[1]} --output-pdf {output}
    '''

##########################################
### Admixture with selection
##########################################

rule admixture_w_ai_bp:
  input:
    expand("{results_dir}/dataframes/admixture.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/admixture_with_incompatability_boxplot.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/admixture_w_ai_bp.py --admixture-pkl {input[0]} --output-pdf {output}
    '''

rule admixture_w_ai_aft:
  input:
    expand("{results_dir}/dataframes/admixture_afts.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/admixture_with_incompatability_aft.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/admixture_aft.py --admixture-afts-pkl {input[0]} --output-pdf {output}
    '''

rule admixture_manhat:
  input:
    expand("{results_dir}/dataframes/admixture.pkl", results_dir=config['results_dir'])
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/admixture_manhattan.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/admixture_manhattan.py --admixture-pkl {input[0]} --output-pdf {output}
    '''

##########################################
### fis
##########################################

rule simulation_fis:
  input:
    expand("{results_dir}/dataframes/fis_pop_split.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/dataframes/fis_mig.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/dataframes/fis_admixture.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/dataframes/fis_tm.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/simulation_fis.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/simulation_fis.py --fis-pop-split-pkl {input[0]} --fis-mig-pkl {input[1]} --fis-admixture-pkl {input[2]} --fis-tm-pkl {input[3]} --output-pdf {output}
    '''

##########################################
### Permutation estimation
##########################################

rule perm_est_fig:
  input:
    expand("{results_dir}/tips_50_perms_100000.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/tips_100_perms_100000.pkl", results_dir=config['results_dir']),
    expand("{results_dir}/tips_200_perms_100000.pkl", results_dir=config['results_dir']),
  resources:
    time="2:30:00",
    threads=1,
    mem="5gb"
  output:
    expand("{figure_dir}/permutation_estimation.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python3 scripts/figures/permutation_estimation.py --tips50 {input[0]} --tips100 {input[1]} --tips200 {input[2]} --output {output}
    '''