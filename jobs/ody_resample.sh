#!/bin/bash

#SBATCH -J c3k_fsps
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 6:00:00 # Runtime
#SBATCH -p conroy_priority,shared,sapphire # Partition to submit to
#SBATCH --mem-per-cpu=2500 # Memory per cpu in MB
#SBATCH -o logs/c3k_resample_%A_%a.log # Standard out goes to this file
#SBATCH -e logs/c3k_resample_%A_%a.log # Standard err goes to this file

# call the script with
# --export=ALL,libname=<libname>,ck_vers=<ck_vers>,synthe=synthe
# with <libname>=(hr/lr/etc.), ck_vers=(c3k_v1.3/c3k_v2.3/etc.), synthe=(vt10_allfal/vt10_optfal etc.)

date; hostname; pwd
module purge
module load python/3.10.9-fasrc01

source activate c3k
PROJECT_DIR=$SCRATCH/conroy_lab/Lab/$USER/c3k-fsps-lib
cd $PROJECT_DIR/src

echo 'libname='$libname
echo 'ck_vers='$ck_vers
echo 'synthe='$synthe

# location of spec and flux HDF5 files
# these should be set when calling the slurm job
#ck_vers=c3k_v2.3
#synthe=vt10_uncal

fulldir="/n/holystore01/LABS/conroy_lab/Lab/bdjohnson/data/kurucz/{}"/${synthe}
# segment specification
segments=$PROJECT_DIR/segments/segments_${libname}.yml

# directory and label for output
seddir=$PROJECT_DIR/output/${ck_vers}/${synthe}/${libname}
mkdir -p $seddir
mkdir -p ${seddir}/for_fsps

zind=${SLURM_ARRAY_TASK_ID}
#zind=-1

python c3k_resample.py --zindex $zind --ck_vers $ck_vers --oldz 0 \
                       --segment_file $segments --oversample 2 --sedname ${libname} \
                       --seddir ${seddir} --fulldir ${fulldir} --specdir spec --fluxdir flux \
                       --bindir ${seddir}/for_fsps \
                       --verbose=False --nowrite 0 --make_seds 1 --make_grid 1 --make_bins 1
