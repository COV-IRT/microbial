#!/bin/bash
#SBATCH -J COVIRT-aws_microbe           # Job name
#SBATCH -o /scratch/06176/jochum00/COVIRT19/stdout/aws.microbe.o%j       # Name of stdout output file
#SBATCH -e /scratch/06176/jochum00/COVIRT19/stderr/aws.microbe.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 64               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 24:00:00        # Run time (hh:mm:ss)
#SBATCH -A COVIRT19       # Allocation name (req'd if you have more than 1)
#SBATCH --qos=vip   #Run as a COVIRT19 VIP in the queue

#this is a reminder on how to pull a docker container using singularity

module load tacc-singularity

singularity run ../../singularity_cache/amazon.sif s3 sync s3://nasa-covid $PWD --exclude "*human*" --endpoint-url=https://s3.wasabisys.com --profile wasabi
