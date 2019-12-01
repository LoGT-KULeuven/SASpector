#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 19:02:15 2019

@author: alerojo
"""

from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import os
import progressbar 

""" kmer

This script calculates the kmers (defined by a provided kmer size) for each missing region, and generates the corresponding abundancy tables and barplots.

"""

def kmer(k, prefix, outdir):
    """ Calculates kmers, generates abundancy table (tsv) and barplot (jpg) which are saved in new subdirectory 'kmer' for every missing region
    
    Parameters
    ----------
    k: int
        k-mer size
    prefix: str
        Name of the genome
    outdir: str
        Output directory
    
    """
    
    newpath = 'kmer'
    os.makedirs(os.path.join(outdir, newpath))
    
    # Define the unmapped regions FASTA file
    bar = progressbar.ProgressBar(widgets = ['Calculating kmers: ', progressbar.Bar(), '(', progressbar.ETA(),')'])
    for i in bar(range(1)):
        unmap = '{outdir}/{prefix}_unmappedregions.fasta'.format(outdir = outdir, prefix = prefix)
    
        for reads in SeqIO.parse(unmap, format = 'fasta'):
        
        # Create kmers and stores them in a dictionary
            kmers = dict()
            for i in range(len(str(reads.seq)) - k+1):
                kmer = str(reads.seq[i:i+k])
                if kmer in kmers:
                    kmers[kmer] += 1
                else:
                    kmers[kmer] = 1
        
        # Write kmers, kmers count and plots
            path = '{outdir}/{k}mer'.format(outdir = outdir, k = str(k))
            with open(os.path.join(path, '{id}_{k}mer.tsv'.format(id = reads.id, k = str(k))), 'w+') as out:
                for key, value in kmers.items():
                    out.write(key + '\t' + str(value) + '\n')
        
        # Create kmer barplots
            plt.figure(figsize = (50, 10))
            sns.set_style('dark')
            fig = sns.barplot(x = list(kmers.keys()), y = list(kmers.values())).set_title('kmer count (k = {k})'.format(k = k))
            #fig = sns.barplot(list(kmers.keys(), list(kmers.values())).set_title('k-mer count (k = {k})'.format(k = k)))
            plt.xlabel('k-mers')
            plt.ylabel('Counts')
            plt.xticks(rotation = 90)
            save = fig.get_figure()
            save.savefig(os.path.join(path, '{id}_{k}mer.jpg'.format(id = reads.id, k = str(k))))
            plt.cla()
            plt.close(save)
            kmers.clear()
