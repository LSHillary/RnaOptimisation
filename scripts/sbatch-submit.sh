#!/bin/bash

#SBATCH --job-name={params.sname}_{rule}
#SBATCH --output=logs/{params.sname}_{rule}%j.out
#SBATCH --error=logs/{params.sname}_{rule}%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={resources.cores}
#SBATCH --mem={resources.mem_mb}MB
#SBATCH --time=10-00:00:00
#SBATCH --partition=high2

srun {exec_job}