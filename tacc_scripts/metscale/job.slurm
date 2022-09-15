#!/bin/bash
#SBATCH -J snakemake           # Job name
#SBATCH -o output.%j       # Name of stdout output file
#SBATCH -e error.%j       # Name of stderr error file
#SBATCH -p skx-normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 04:00:00        # Run time (hh:mm:ss)
#SBATCH -A COVIRT19       # Allocation name (req'd if you have more than 1)
#SBATCH --qos=vip

#SBATCH --mail-user=kternus@signaturescience.com
#SBATCH --mail-type=all    # Send email at begin and end of job

# load modules first
module load tacc-singularity
module list

# set umask to enable sharing across the COVIRT19 group
umask 0007

conda activate metag

export SINGULARITY_BINDPATH="data:/tmp"
snakemake --cores --use-singularity --configfile=config/my_custom_config.json tax_class_bracken_workflow > snakemake.${SLURM_JOBID}.log 2>&1

