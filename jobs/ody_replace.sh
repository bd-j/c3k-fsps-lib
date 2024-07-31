#!/bin/bash

#SBATCH -J c3k_replacefal
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1:00:00 # Runtime
#SBATCH -p conroy_priority,shared # Partition to submit to
#SBATCH --mem-per-cpu=3GB # Memory per node in MB (see also --mem-per-cpu)
#SBATCH --output logs/c3k_replace_%A_%a.log # Standard out goes to this file

date; hostname; pwd

module purge
module load python/3.10.9-fasrc01

source activate c3k
PROJECT_DIR=$SCRATCH/conroy_lab/$USER/c3k-fsps-lib
cd $PROJECT_DIR/src

# location of spec and flux HDF5 files
ck_vers=c3k_v2.3
synthe=vt10_uncal
#fulldir="/n/holystore01/LABS/conroy_lab/Lab/bdjohnson/data/kurucz/{}"/${synthe}


python replace_optfal.py --zindex ${SLURM_ARRAY_TASK_ID} --synthe $synthe --ck_vers $ck_vers
