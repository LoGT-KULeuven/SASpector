#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 16:40:09 2019

@author: alerojo
"""

import pandas as pd
import pylab
from Bio import SeqIO
from Bio.SeqUtils import GC

# Mauve - Parsing positions
backbone = pd.read_table("positions.backbone", sep= "\s+")
regions = backbone[(backbone.seq1_leftend == 0) & (backbone.seq1_rightend == 0)]
regions = regions[["seq0_leftend","seq0_rightend"]]            

missing = {}    
for records in SeqIO.parse("assembly.fasta","fasta"):
    for i in range(0, regions.shape[0]):
        start = regions.iloc[i,0]
        end = regions.iloc[i,1]
        missing[i] = str(records.seq[start:end])

# Fasta file - Missing regions
with open("MissingRegions.fasta","w+") as fasta:
    for keys in missing:
        fasta.write(">"+str(keys)+"\n"+missing[keys]+"\n")
        
# Sequence Analysis
gc_values = sorted(GC(rec.seq) for rec in SeqIO.parse("MissingRegions.fasta","fasta"))
sizes = [len(rec) for rec in SeqIO.parse("MissingRegions.fasta","fasta")]

# GC - Plot
pylab.plot(gc_values)
pylab.title("%i Missing regions sequences\nGC%% %0.1f to %0.1f" \
            % (len(gc_values), min(gc_values), max(gc_values)))
pylab.xlabel("Regions")
pylab.ylabel("GC%")
pylab.show()

# Length - Plot
pylab.hist(sizes, bins=20)
pylab.title("%i Missing regions sequences\nLengths %i to %i" \
            % (len(sizes), min(sizes), max(sizes)))
pylab.xlabel("Sequence length (bp)")
pylab.ylabel("Count")
pylab.show()