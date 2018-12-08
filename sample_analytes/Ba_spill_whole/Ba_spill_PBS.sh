#!/bin/bash

#PBS -A open
#PBS -l nodes=1:ppn=20
#PBS -l walltime=48:00:00
#PBS -l pmem=10gb
#PBS -N Ba_spill

cd $PBS_O_WORKDIR

module load r/3.4

Rscript Ba_spill_clust.R
