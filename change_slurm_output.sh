#!/bin/bash

VARS_PER_GENE=(2 3 5 7 10)

READS_PER_SITE=(10 20 50 100)

THETA=(0.25 0.5 0.75)

i=1

for var in ${VARS_PER_GENE[@]}; do
    for reads in ${READS_PER_SITE[@]}; do
        for thet in ${THETA[@]}; do
            cp /datacommons/allenlab/stephanie/sim-slurms/outputs/${i}.output /datacommons/allenlab/stephanie/trio_output/var_${var}_reads_${reads}_theta_${thet}.output
	    i=$((i+1))
	done
    done
done
