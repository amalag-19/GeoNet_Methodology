#!/bin/bash

#PBS -A open
#PBS -l nodes=1:ppn=5
#PBS -l walltime=24:00:00
#PBS -l pmem=2gb

cd $PBS_O_WORKDIR

module load r/3.4

#LD_PRELOAD="libmkl_gf_lp64.so libmkl_sequential.so libmkl_core.so" R --no-save --slave --args $args < runningcode.R

R --no-save --slave --args $args < clrun_Br_projected.R
