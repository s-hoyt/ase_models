#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2021 William H. Majoros <bmajoros@alumni.duke.edu>
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import os
import sys
from Pipe import Pipe
import ProgramName
from Rex import Rex
rex=Rex()

def parseParm(parmString,label):
    rex.findOrDie(label+"_([\d\.]+)",parmString)
    return rex[1]

def parseParmString(parmString):
    sites=int(parseParm(parmString,"var"))
    reads=int(parseParm(parmString,"reads"))
    theta=float(parseParm(parmString,"theta"))
    return (sites,reads,theta)

def evalTrioModel(parmString,dataDir,outputsDir):
    truth=dataDir+"/"+parmString+"_truth.essex"
    predictions=outputsDir+"/"+parmString+".output"
    cmd="python3 ~/ase_models/eval1_stephanie.py "+truth+" "+predictions
    text=Pipe.run(cmd)
    rex.findOrDie("(.*)% correct among trios with evidence",text)
    acc=float(rex[1])/100
    return acc
    
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <data-dir> <outputs-dir>\n")
(dataDir,outputsDir)=sys.argv[1:]

print("sites\treads\ttheta\tacc")
files=os.listdir(dataDir)
for filename in files:
    if(not rex.find("(.+)_truth.essex",filename)): continue
    parmString=rex[1]
    #print(filename)
    #print(parmString)
    (sites,reads,theta)=parseParmString(parmString)
    #print(sites,reads,theta)
    trioAcc=evalTrioModel(parmString,dataDir,outputsDir)
    print(sites,reads,theta,round(trioAcc,2),sep="\t")
    

#var_2_reads_1000_theta_0.75_truth.essex
