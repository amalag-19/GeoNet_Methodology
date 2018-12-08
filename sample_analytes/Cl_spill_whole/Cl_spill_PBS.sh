#!/bin/bash

#PBS -A open
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -l pmem=10gb
#PBS -N Cl_spill

cd $PBS_O_WORKDIR

module load r/3.4

Rscript Cl_spill_clust.R
