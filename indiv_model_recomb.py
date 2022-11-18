#!/usr/bin/env python

import os

code_dir="/hpc/home/sh700/TrioBEAST/"
model_dir="/hpc/home/sh700/ase_models/"
data_dir="/hpc/group/majoroslab/stephanie/"

#RECOMB=(0, 0.0001, 0.001, 0.005, 0.01, 0.1, 0.5, 1)
RECOMB=[.001]
MCMC = 1000
LAMBDA = 1.2


for r in RECOMB:
    for indiv in ["child", "mother", "father"]:
        py = "python " + code_dir + "run-phased-indiv.py " + model_dir + "PhasedIndividual "
        params = data_dir + "data_sim_thresh/recomb_" + str(r) + "_phased.essex " + indiv + " " + str(MCMC) + " " + str(LAMBDA)
        output = data_dir + "indiv_output_thresh/" + indiv + "_recomb_" + str(r) + ".txt"
        cmd = py + params + " > " + output
        os.system(cmd)
