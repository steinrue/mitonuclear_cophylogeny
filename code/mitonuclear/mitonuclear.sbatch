#!/bin/bash
#SBATCH --job-name=mitonuclear
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=240:00:00
#SBATCH --mem=1gb
#SBATCH --output=/gpfs/data/steinruecken-lab/users/rowan/github_test/slurm_tmp/%x.o%j
#SBATCH --error=/gpfs/data/steinruecken-lab/users/rowan/github_test/slurm_tmp/%x.e%j
snakemake --rerun-incomplete --keep-going --cluster 'sbatch --nodes=1 --ntasks=1 --cpus-per-task={resources.threads} --mem={resources.mem} --time={resources.time} --job-name={rule} --output=/gpfs/data/steinruecken-lab/users/rowan/github_test/snakeout_tmp/%j/%x.o%j --error=/gpfs/data/steinruecken-lab/users/rowan/github_test/snakeout_tmp/%j/%x.e%j' --jobs 2000