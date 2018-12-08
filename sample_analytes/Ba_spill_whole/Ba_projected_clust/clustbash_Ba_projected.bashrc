#!/bin/bash
for args in `seq 1 80`; 
do
	qsub clustshell_Ba_projected.sh -v args=$args
done
