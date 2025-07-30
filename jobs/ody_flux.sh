#!/bin/bash

#SBATCH -J c3k_fluxfiles
#SBATCH -n 20 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 6:00:00 # Runtime
#SBATCH -p conroy # Partition to submit to
#SBATCH --mem-per-cpu=2000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH -o logs/c3k_flux_%A_%a.out # Standard out goes to this file
#SBATCH -e logs/c3k_flux_%A_%a.err # Standard err goes to this file


source activate c3k
mkdir -p $SCRATCH/conroy_lab/Lab/$USER/run_ckc
cd $SCRATCH/conroy_lab/$USER/c3k-fsps-lib/src

python make_flux.py --np=${SLURM_JOB_CPUS_PER_NODE} --feh=-99 --ck_vers=c3k_v1.3