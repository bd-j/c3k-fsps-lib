#!/bin/bash

#SBATCH -J c3k_copy
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1:00:00 # Runtime
#SBATCH -p conroy_priority,shared # Partition to submit to
#SBATCH --mem-per-cpu=4000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH -o logs/c3k_copy_%A.log # Standard out goes to this file
#SBATCH -e logs/c3k_copy_%A.log # Standard err goes to this file

date; hostname; pwd
echo 'libname='$libname


# directory and label for output
PROJECT_DIR=$SCRATCH/conroy_lab/$USER/c3k-fsps-lib
cvers=c3k_v2.3
synthe=vt10_uncal
seddir=$PROJECT_DIR/output/${cvers}_${synthe}/${libname}
storage=/n/holystore01/LABS/conroy_lab/Lab/bdjohnson/data/kurucz/${cvers}/${synthe}/fsps-lib
mkdir -p $storage
rm -rf ${storage}/${libname}
cp -r ${seddir} ${storage}/
