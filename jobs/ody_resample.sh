#!/bin/bash

#SBATCH -J c3k_fsps
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 6:00:00 # Runtime
#SBATCH -p conroy_priority,shared # Partition to submit to
#SBATCH --mem-per-cpu=2500 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH -o logs/c3k_resample_%A_%a.log # Standard out goes to this file
#SBATCH -e logs/c3k_resample_%A_%a.log # Standard err goes to this file

# call the script with --export=ALL,libname=<libname> (hr/lr/etc.)

date; hostname; pwd
echo $libname

module purge
module load python/3.10.9-fasrc01

source activate c3k
PROJECT_DIR=$SCRATCH/conroy_lab/$USER/c3k-fsps-lib
cd $PROJECT_DIR/src

# location of spec and flux HDF5 files
ck_vers=c3k_v2.3
synthe=optfal
fulldir="/n/holystore01/LABS/conroy_lab/Lab/bdjohnson/data/kurucz/{}"/${synthe}
# segment specification
segments=$PROJECT_DIR/segments/segments_${libname}.yml

# directory and label for output
seddir=$PROJECT_DIR/output/${ck_vers}_${synthe}/${libname}
mkdir -p $seddir
mkdir -p ${seddir}/for_fsps

python c3k_resample.py --zindex ${SLURM_ARRAY_TASK_ID} --ck_vers $ck_vers --oldz 0 \
                       --segment_file $segments --oversample 2 --sedname ${libname} \
                       --seddir ${seddir} --fulldir ${fulldir} --specdir spec --fluxdir flux \
                       --bindir ${seddir}/for_fsps \
                       --verbose=False --nowrite 0 --make_seds 1 --make_grid 1 --make_bins 1
