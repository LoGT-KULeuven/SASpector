#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 11:46:33 2019

@author: alerojo
"""       

import subprocess
import shlex
import progressbar 
#from Bio.Blast.Applications import NcbiblastxCommandline as blastx

def prokka(prefix, outdir):
    bar = progressbar.ProgressBar(widgets = ['Predicting genes: ', progressbar.Bar(), '(', progressbar.ETA(),')'])
    for i in bar(range(1)):
        cmd = 'prokka --outdir {outdir}/genesprediction --prefix {prefix}.predictedgenes {outdir}/{prefix}_unmappedregions.fasta'.format(prefix = prefix, outdir = outdir)
        process = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
        while process.poll() is None:
            l = process.stdout.readline()
    
    #blastx -query Aba1_predictedgenes.fsa -subject saspector_proteindb.fasta -evalue 0.001 -outfmt '6 qseqid qstart qend sseqid sstart send pident evalue qcovs'
#def blast():
#    cline = blastx(cmd = 'blastx', out = '{prefix}_blastxresults.tsv'.format(prefix = prefix), evalue = 0.001, 
#                   outfmt = '6 qseqid qstart qend sseqid sstartstdout, stderr = cline()stdout, stderr = cline() send pident evalue qcovs', query = )
#    stdout, stderr = cline()
 

    
