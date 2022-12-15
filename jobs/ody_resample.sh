#!/bin/bash

#SBATCH -J c3k_fsps
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 6:00:00 # Runtime
#SBATCH -p conroy # Partition to submit to
#SBATCH --mem-per-cpu=2000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH -o /n/holyscratch01/conroy_lab/bdjohnson/alpha-fsps/jobs/logs/c3k_R10k_%A_%a.out # Standard out goes to this file
#SBATCH -e /n/holyscratch01/conroy_lab/bdjohnson/alpha-fsps/jobs/logs/c3k_R10k_%A_%a.err # Standard err goes to this file

module purge
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02 hdf5/1.10.1-fasrc01
module load python/2.7.14-fasrc01

source activate pro
PROJECT_DIR=$SCRATCH/conroy_lab/$USER/c3k-fsps-lib
cd $PROJECT_DIR/src

# location of fullres and flux HDF5 files
fulldir=/n/holystore01/LABS/conroy_lab/Lab/bdjohnson/data/kurucz/{}/fullres
# segment specification
segments=$PROJECT_DIR/segements/segments_${libname}.yml
# directory and label for output
seddir=$PROJECT_DIR/output/${libname}
mkdir -p $seddir

python c3k_resample.py --zindex ${SLURM_ARRAY_TASK_ID} --ck_vers c3k_v1.3 \
                       --segment_file $segments --oversample 2 \
                       --seddir ${seddir} --sedname ${libname} --fulldir ${fulldir} \
                       --verbose=False --nowrite 1