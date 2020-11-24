#!/bin/bash
#SBATCH -J kneaddata           # Job name
#SBATCH -o /scratch/06176/jochum00/COVIRT19/stdout/knead.o%j       # Name of stdout output file
#SBATCH -e /scratch/06176/jochum00/COVIRT19/stderr/knead.e%j       # Name of stderr error file
#SBATCH -p normal       # Queue (partition) name
#SBATCH -N 80               # Total # of nodes (must be 1 for serial)
#SBATCH -n 80               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 24:00:00        # Run time (hh:mm:ss)
#SBATCH -A COVIRT19       # Allocation name (req'd if you have more than 1)
#SBATCH --qos=vip         #Run as a COVIRT19 VIP in the queue



#Load the modules you need here
module load launcher
module load tacc-singularity

umask 0007
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=commands.txt
export LAUNCHER_SCHED=dynamic
export LAUNCHER_PPN=1

$LAUNCHER_DIR/paramrun

exit 00