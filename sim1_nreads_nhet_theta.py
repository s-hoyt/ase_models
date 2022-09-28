#!/usr/bin/env python

#=========================================================================

# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public

# License (GPL) version 3, as described at www.opensource.org.

# Copyright (C)2017 William H. Majoros (martiandna@gmail.com).

#=========================================================================

from __future__ import (absolute_import, division, print_function,

   unicode_literals, generators, nested_scopes, with_statement)

from builtins import (bytes, dict, int, list, object, range, str, ascii,

   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)

# The above imports should allow this program to run in both Python 2 and

# Python 3.  You might need to update your version of module "future".

import sys

from SlurmWriter import SlurmWriter

#arguments to code/sim1 in order

VCF="/datacommons/allenlab/ThousandGenomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

PARENT_1="HG00097"
PARENT_2="HG00096"

#make this variable for now to test python script,
#but this is probably not an interseting paramter to change for actual analysis
NGENES=1000

VARS_PER_GENE=(2,3,5,7,10)

READS_PER_SITE=(10,20,50,1000)

RECOMB=0.001

THETA=(0.25,0.5,0.75)

MCMC=1000

#vars to build base of command

BASE="/datacommons/allenlab/stephanie/"

SLURM_DIR=BASE+"sim-slurms/"

OUT_DIR=BASE+"data_sim/"

JOB_NAME="SIM"

MAX_PARALLEL=1000

#ADDITIONAL_LINES="#SBATCH --exclude=x1-01-1,x1-01-2"

ADDITIONAL_LINES=""

def makeCommand(vcf,p1,p2,ngenes,nvars,reads,recomb,theta,mcmc):

    truth = OUT_DIR + "var_" + str(nvars) + "_reads_" + str(reads) + "_theta_" + str(theta) + "_truth.essex"
    data = OUT_DIR + "var_" + str(nvars) + "_reads_" + str(reads) + "_theta_" + str(theta) + "_data.essex"

    sim_cmd = "code/sim1" + " " + vcf + " " + p1 + " " + p2 + " " + str(ngenes) + \
              " " + str(nvars) + " " + str(reads) + " " + \
              str(recomb) + " " + str(theta) + " " + truth + " " + data

    phased = OUT_DIR + "var_" + str(nvars) + "_reads_" + str(reads) + "_theta_" + str(theta) + "_phased.essex"

    phase_cmd = "code/phase-trio" + " " + data + " " + phased

    trio_cmd = "LAMBDA=1.2; AFFECTED=0.5; " +\
                "code/refactored.py -s stan.out code/Refactored " + \
                phased + " " + str(mcmc) + " " + \
                "0-" + str(ngenes-1) + " " + "$LAMBDA $AFFECTED"

    cmd = "cd " + BASE + "\n\n" + \
           sim_cmd + "\n\n" + \
           phase_cmd + "\n\n" + \
           trio_cmd

   
    return cmd



def addCommands(slurm, var, reads, thet):

    slurm.addCommand(makeCommand(VCF, PARENT_1, PARENT_2,
                                 NGENES, var, reads,
                                 RECOMB, thet, MCMC))


#=============================================================

# main()

#=============================================================

slurm=SlurmWriter()

for var in VARS_PER_GENE:
    for reads in READS_PER_SITE:
        for thet in THETA:
            addCommands(slurm, var, reads, thet)

slurm.mem(1500)

slurm.setQueue("scavenger")

slurm.writeArrayScript(SLURM_DIR,JOB_NAME,MAX_PARALLEL,ADDITIONAL_LINES)
