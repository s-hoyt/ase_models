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
import sys
import os
import ProgramName
from EssexParser import EssexParser
from Rex import Rex
rex=Rex()

class BufferedReader:
    def __init__(self,filename):
        self.fh=open(filename,"rt")
        self.fh.readline() # discard the header line
        self.buffer=None
    def nextPred(self):
        if(self.buffer is not None):
            temp=self.buffer
            self.buffer=None
            return temp
        while(True):
            line=self.fh.readline()
            while(not rex.find("^GENE",line)): 
                line=self.fh.readline()
                if(line is None): return None
            print(line)
            fields=line.rstrip().split()
            #(ID,Palt,theta,CI,Pchild)=fields
            (ID,Palt,theta,CI)=fields
            line=self.fh.readline()
            if(not rex.find("(\d)% : (\d)(\d) (\d)(\d) (\d)(\d)",line)):
                raise Exception("Can't parse: "+line)
            posterior=rex[1]; mother=[rex[2],rex[3]];
            father=[rex[4],rex[5]]; child=[rex[6],rex[7]]
            return [ID,Palt,theta,posterior,mother,father,child]
          
def countHets(root,indiv):
    sites=root.findChildren("site")
    num_het=[]
    for site in sites:
        node=site.findChild("genotypes").findChild(indiv)
        if(node is None): raise Exception("Can't find node in hasHets")
        num_het+=[int(node[0]!=node[1])]
    return num_het

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <truth.essex dir> \n")

this_dir = sys.argv[1:][0]
all_data_files = os.listdir(this_dir)
truth_files = [x for x in all_data_files if "truth" in x]


for truthFile in truth_files:

    essexReader=EssexParser(this_dir + truthFile)
    numGenes_trihet=0;
    numSites_trihet=0
    while(True):
        root=essexReader.nextElem()
        if(root is None): break
    
        res_m = countHets(root,"mother")
        res_f = countHets(root,"father")
        res_c = countHets(root,"child")
    
        res_trip = [res_m[i] + res_f[i] + res_c[i] for i in range(len(res_m))]
        num_tri_hets = res_trip.count(3)

        numSites_trihet += num_tri_hets
        numGenes_trihet += int(num_tri_hets >= 1)

    print(truthFile)
    print(numGenes_trihet," number of genes with at least 1 triple het site out of 1000 genes \n")
    print(numSites_trihet," number of triple het sites\n")
