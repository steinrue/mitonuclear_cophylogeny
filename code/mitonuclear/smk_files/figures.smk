configfile: "config.yaml"

rule one_kg_fis:
  input:
    expand("{results_dir}/zscores/whole_genome.zscores", results_dir=config['results_dir']),
    expand("{results_dir}/FIS/all_populations.fis", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/1kg_fis.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/1kg_fis.py --whole_genome_zscores {input[0]} --all_populations_fis {input[1]} --output_figure {output}
    '''

rule one_kg_cophylogeny_scores:
  input:
    expand("{results_dir}/zscores/whole_genome.zscores", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/1kg_cophylogeny_scores.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/1kg_cophylogeny_scores.py --whole_genome_zscores {input[0]} --output_figure {output}
    '''
  
rule chr20_vs_fis:
  input:
    expand("{results_dir}/zscores/chr20.zscores", results_dir=config['results_dir']),
    expand("{results_dir}/FIS/all_populations.fis", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/1kg_chr20_vs_fis.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/1kg_chr20_vs_fis.py --chr20_zscores {input[0]} --all_populations_fis {input[1]} --output_figure {output}
    '''

rule single_population_enrichment:
  input:
    expand("{results_dir}/single_locus/enrichment_df.pkl", results_dir=config['results_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/enrichment.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/enrichment.py --enrichment_df {input[0]} --output_pdf {output}
    '''

rule get_manhattan:
  input:
    expand("{results_dir}/single_locus/{POP}.pvals", results_dir=config['results_dir'], POP = "{POP}"),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/{POP}_manhattan.pdf", figure_dir=config['figure_dir'], POP = "{POP}")
  shell:
    '''
    python scripts/figures/get_manhattan.py --single_locus_pvals {input[0]} --output_figure {output}
    '''

rule get_peaks:
  input:
    expand("{results_dir}/zscores/whole_genome.zscores", results_dir=config['results_dir']),
    expand("{data_dir}/MitoCarta/Human.MitoCarta3.0.xls", data_dir=config['data_dir']),
    expand("{data_dir}/GenCode/gencode.v49.annotation.gff3.gz", data_dir=config['data_dir']),
    expand("{results_dir}/single_locus/JPT.pvals", results_dir=config['results_dir'], POP = "{POP}"),
    expand("{results_dir}/single_locus/IBS.pvals", results_dir=config['results_dir'], POP = "{POP}"),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/peaks.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/get_peaks.py --whole_genome_zscores {input[0]} --gencode_gff3_gz {input[2]} --mitocarta_xls {input[1]} --jpt_pvals {input[3]} --ibs_pvals {input[4]} --output_pdf {output}
    '''

rule get_jpt_mt_tree:
  input:
    expand("{data_dir}/GenCode/gencode.v49.annotation.gff3.gz", data_dir=config['data_dir']),
    expand("{results_dir}/single_locus/JPT.pvals", results_dir=config['results_dir']),
    expand("{tmp_dir}/JPT_chr19.vcf.gz", tmp_dir=config['tmp_dir']),
    expand("{tmp_dir}/JPT_MT.hsd", tmp_dir=config['tmp_dir']),
    expand("{trees_dir}/JPT/JPT_chr19.trees", trees_dir=config['trees_dir']),
    expand("{tmp_dir}/JPT_samples.list", tmp_dir=config['tmp_dir']),
    expand("{trees_dir}/JPT/JPT_chrMT.nwk", trees_dir=config['trees_dir']),
  resources:
    time="0:30:00",
    threads=1,
    mem="1gb"
  output:
    expand("{figure_dir}/JPT_MT_tree.pdf", figure_dir=config['figure_dir'])
  shell:
    '''
    python scripts/figures/jpt_mt_tree.py --gencode_path {input[0]} --results_df {input[1]} --vcf {input[2]} --hsd {input[3]} --trees {input[4]} --samples_list {input[5]} --mt_tree {input[6]} --output_pdf {output}
    '''