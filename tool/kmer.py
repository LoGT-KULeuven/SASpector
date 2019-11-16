#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 19:02:15 2019

@author: alerojo
"""

from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess
import shlex
import os

newpath = 'kmer'
os.makedirs(newpath)
global newpath

def kat(k, reference, prefix):
    cmd = 'kat hist {reference} -o {prefix} -m {k}'.format(reference = reference, prefix = prefix, k = k)
    process = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE)
    while process.poll() is None:
        l = process.stdout.readline()
        print(l.decode())
    print(process.stdout.read().decode())

def kmer(k, prefix):
    
    # Define the unmapped regions FASTA file
    unmap = '{prefix}_unmappedregions.fasta'
    
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
        with open(os.path.join(newpath, '{id}_kmer.tsv'.format(id = reads.id)), 'w+') as out:
            for key, value in kmers.items():
                out.write(key + '\t' + str(value), + '\n')
        
        # Create kmer barplots
        plt.figure(figsize = (50, 10))
        sns.set_style('dark')
        fig = sns.barplot(list(kmers.keys(), list(kmers.values())).set_title('k-mer count (k = {k})'.format(k = k)))
        plt.xlabel('k-mers')
        plt.ylabel('Counts')
        plt.xticks(rotation = 90)
        fig.savefig(os.path.join(newpath, '{id}_kmer.jpg'.format(id = reads.id)))
        kmers.clear()
    
        




    
