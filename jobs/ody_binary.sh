#!/bin/bash

#SBATCH -J c3k_fs_bin
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 00:30:00 # Runtime
#SBATCH -p conroy # Partition to submit to
#SBATCH --mem-per-cpu=2000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH -o logs/c3k_bin_%A_%a.log # Standard out goes to this file
#SBATCH -e logs/c3k_bin_%A_%a.log # Standard err goes to this file

module purge
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02 hdf5/1.10.1-fasrc01
module load python/2.7.14-fasrc01

source activate c3k
PROJECT_DIR=$SCRATCH/conroy_lab/$USER/c3k-fsps-lib
cd $PROJECT_DIR/src

# directory with fsps files
seddir=$PROJECT_DIR/output/${libname}
#directory for the output
outdir=${seddir}/for_fsps
mkdir -p $outdir

# Note there are 5 afes
#ii=( 0 1 2 3 4 )
# solar is ii=1

python c3k_binary.py --zindex=${SLURM_ARRAY_TASK_ID} --ck_vers=c3k_v1.3 \
                     --seddir=${seddir} --sedname=${libname} --outdir=${outdir}

wc ${outdir}/*.lambda
wc ${outdir}/*.zlegend.dat