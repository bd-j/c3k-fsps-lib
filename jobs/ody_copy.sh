#!/bin/bash

#SBATCH -J c3k_copy
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1:00:00 # Runtime
#SBATCH -p conroy_priority,shared,sapphire # Partition to submit to
#SBATCH --mem-per-cpu=4000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH -o logs/c3k_copy_%A.log # Standard out goes to this file
#SBATCH -e logs/c3k_copy_%A.log # Standard err goes to this file

# call the script with
# --export=ALL,libname=<libname>,ck_vers=<ck_vers>,synthe=synthe
# with <libname>=(hr/lr/etc.), ck_vers=(c3k_v1.3/c3k_v2.3/etc.), synthe=(vt10_allfal/vt10_optfal etc.)

date; hostname; pwd

echo 'libname='$libname
echo 'ck_vers='$ck_vers
echo 'synthe='$synthe

# directory and label for output
PROJECT_DIR=$SCRATCH/conroy_lab/$USER/c3k-fsps-lib
seddir=$PROJECT_DIR/output/${ck_vers}/${synthe}/${libname}
storage=/n/holystore01/LABS/conroy_lab/Lab/bdjohnson/data/kurucz/${ck_vers}/${synthe}/fsps-lib

mkdir -p $storage
rm -rf ${storage}/${libname}
cp -r ${seddir} ${storage}/
