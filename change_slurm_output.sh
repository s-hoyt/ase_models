#!/bin/bash

VARS_PER_GENE=(2 3 5 7)

READS_PER_SITE=(10 20 50 100)

THETA=(0.25 0.75)

i=1

for var in ${VARS_PER_GENE[@]}; do
    cp /hpc/group/majoroslab/stephanie/sim_triplehets/outputs/${i}.output /hpc/group/majoroslab/stephanie/no_triphet_trio_output/var_${var}_reads_30_theta_0.5.output
    i=$((i+1))
done
    
for reads in ${READS_PER_SITE[@]}; do
    cp /hpc/group/majoroslab/stephanie/sim_triplehets/outputs/${i}.output /hpc/group/majoroslab/stephanie/no_triphet_trio_output/var_10_reads_${reads}_theta_0.5.output
    i=$((i+1))
done

    
for thet in ${THETA[@]}; do
    cp /hpc/group/majoroslab/stephanie/sim_triplehets/outputs/${i}.output /hpc/group/majoroslab/stephanie/no_triphet_trio_output/var_10_reads_30_theta_${thet}.output
    i=$((i+1))
done
