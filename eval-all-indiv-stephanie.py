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
    #sites=int(parseParm(parmString,"sites"))
    #reads=int(parseParm(parmString,"reads"))
    #theta=float(parseParm(parmString,"theta"))
    recomb=float(parseParm(parmString, "recomb"))
    return (recomb)

def evalModel(parmString,dataDir,outputsDir,threshold):
    truth=dataDir + parmString + "_truth.essex"
    cmd="python ~/ase_models/eval-indiv-stephanie.py "+truth+" "+outputsDir+" "+threshold
    text=Pipe.run(cmd)
    rex.findOrDie("(.*)% correct among trios with evidence",text)
    acc=float(rex[1])/100.0
    return acc
    
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <data-dir> <outputs-dir> <p-threshold>\n")
(dataDir,outputsDir,threshold)=sys.argv[1:]

#print("sites\treads\ttheta\tacc")
print("recomb\tacc")
files=os.listdir(dataDir)
for filename in files:
    if(not rex.find("(.+)_truth.essex",filename)): continue
    parmString=rex[1]
    (recomb)=parseParmString(parmString)
    trioAcc=evalModel(parmString,dataDir,outputsDir,threshold)
    print(recomb,round(trioAcc,2),sep="\t")
    
