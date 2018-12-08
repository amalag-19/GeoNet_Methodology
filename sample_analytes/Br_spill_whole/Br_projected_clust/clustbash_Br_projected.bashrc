#!/bin/bash
for args in `seq 1 100`; 
do
	qsub clustshell_Br_projected.sh -v args=$args
done
