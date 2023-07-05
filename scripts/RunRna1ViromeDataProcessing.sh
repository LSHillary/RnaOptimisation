#!/bin/bash
#SBATCH --job-name=RNA1snakemake
#SBATCH --output=logs/RNA1snakemake_%j.out
#SBATCH --error=logs/RNA1snakemake_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=lhillary@ucdavis.edu
#SBATCH --job-name=sortmerna
#SBATCH --nodes=1
#SBATCH -t 10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=high2
source ~/.bashrc
cd rna1
micromamba activate ViromeDataProcessing
snakemake --snakefile ../scripts/0-preprocessing.smk --profile slurm --configfile ../rna1_pipeline_config.yml
#snakemake --snakefile ../scripts/1.1-QC.smk --profile slurm
#snakemake --snakefile ../scripts/sortmerna.smk --profile slurm
#snakemake --snakefile ../scripts/1.2-bbcms.smk --profile slurm
#snakemake --snakefile ../scripts/1.2-BBNorm.smk --profile slurm
#snakemake --snakefile ../scripts/megahit.smk --profile slurm
