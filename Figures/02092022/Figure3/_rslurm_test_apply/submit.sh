#!/bin/bash
#
#SBATCH --array=0-1
#SBATCH --cpus-per-task=2
#SBATCH --job-name=test_apply
#SBATCH --output=slurm_%a.out
C:/PROGRA~1/R/R-41~1.2/bin/x64/Rscript --vanilla slurm_run.R
