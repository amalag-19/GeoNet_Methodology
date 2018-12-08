#!/bin/bash
for args in `seq 1 80`; 
do
	qsub clustshell_Mg_projected.sh -v args=$args
done
