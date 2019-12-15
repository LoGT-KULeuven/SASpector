#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 11:46:33 2019

@author: alerojo
"""       

import subprocess
import shlex
import progressbar 
from Bio.Blast.Applications import NcbiblastxCommandline as blastx

""" gene_predict

This script predicts genes that are in the missing regions of the reference FASTA file using Prokka and performs
a blastx to search hits across several bacterial species from an in-house protein database using BLAST+. The input file for Prokka
is the missing regions FASTA file generated in the summary script. The input files for blastx are the predicted genes DNA sequences FASTA file
and the bacterial protein sequences in FASTA file.

"""

def prokka(prefix, outdir):
    """ Wraps Prokka and generates the predicted genes per missing region with their annotation
    
    Parameters
    ----------
    prefix: str
        Name of the genome
    outdir: str
        Output directory
        
    """
    bar = progressbar.ProgressBar(widgets = ['Predicting genes: ', progressbar.Bar(), '(', progressbar.ETA(),')'])
    for i in bar(range(1)):
        cmd = 'prokka --outdir {outdir}/genesprediction --prefix {prefix}.predictedgenes {outdir}/{prefix}_unmappedregions.fasta'.format(prefix = prefix, outdir = outdir)
        process = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
        while process.poll() is None:
            l = process.stdout.readline() 
            
def blast(outdir, prefix, proteindb):
    """ Wraps BLAST+ and performs blastx with output Prokka nucleotide FASTA file and SASpector protein FASTA file
    
    Parameters
    ----------
    prefix: str
        Name of the genome
    outdir: str
        Output directory
    proteindb: str
        The file location of the protein database FASTA file
    
    """
    bar = progressbar.ProgressBar(widgets = ['BLAST genes: ', progressbar.Bar(), '(', progressbar.ETA(),')'])
    for i in bar(range(1)):
        cline = blastx(cmd = 'blastx', out = '{outdir}/{prefix}_blastxresults.tsv'.format(prefix = prefix, outdir = outdir), evalue = 0.001, 
                   outfmt = '6 qseqid qstart qend sseqid sstartstdout, stderr = cline()stdout, stderr = cline() send pident evalue qcovs', query = '{outdir}/genesprediction/{prefix}.predictedgenes.fsa'.format(outdir = outdir, prefix = prefix), subject = proteindb)
        stdout, stderr = cline()

    
