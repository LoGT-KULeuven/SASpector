#!python
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
import shlex
import shutil
import subprocess
import glob
import sourmash
import pandas as pd
import numpy as np


""" kmer

This script calculates the kmers for each missing regions defined by a kmer size and generate abundances table and barplot
for the kmers. Additionally, it runs Tandem Repeat Finder for tandem repeats and creates pairwisie comparison using sourmash.

"""

def kmer(k, prefix, outdir):
    """ Wraps Prokka and generates the predicted genes per missing region with their annotation
    
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
    os.makedirs(os.path.join(outdir,newpath))
    
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
            path = '{outdir}/kmer'.format(outdir = outdir)
            with open(os.path.join(path, '{id}_kmer.tsv'.format(id = reads.id)), 'w+') as out:
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
            save.savefig(os.path.join(path, '{id}_kmer.jpg'.format(id = reads.id)))
            plt.cla()
            plt.close(save)
            kmers.clear()
    
def trf(prefix,outdir):
    """ Wraps TRF and generates tandem repeats reports for the missing regions
    
    Parameters
    ----------
    prefix: str
        Name of the genome
    outdir: str
        Output directory
    
    """
    
    cmd = 'trf409.linux64 {outdir}/{prefix}_unmappedregions.fasta 2 5 7 80 10 50 2000'.format(prefix = prefix, outdir = outdir)
    process = subprocess.run(shlex.split(cmd), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    newdir = 'trf'
    os.makedirs(os.path.join(outdir, newdir))
    dest = os.path.join(outdir,newdir)
    for html in glob.glob('*.html'):
        shutil.move(html,dest)
        
        
def clustermap(prefix, outdir):
    """ Computes the pairwise comparison between kmers (k = 31) of missing regions and mapped regions using Jaccard similarity.
    Finally, generates a cluster map for those comparisons.
    
    Parameters
    ----------
    prefix: str
        Name of the genome
    outdir: str
        Output directory

    """
    
    regions_fasta = ['{outdir}/{prefix}_unmappedregions.fasta'.format(outdir = outdir, prefix = prefix),
                     '{outdir}/{prefix}_mappedregions.fasta'.format(outdir = outdir, prefix = prefix)]
    minhashes = list()
    id_records = list()
    
    for r in regions_fasta:
        E = sourmash.MinHash(n = 1000, ksize = 31)
        for record in SeqIO.parse(r, format = 'fasta'):
            E.add_sequence(str(record.seq))
            if r == '{outdir}/{prefix}_unmappedregions.fasta'.format(outdir = outdir, prefix = prefix):
                newid = ''.join([record.id, '_Um'])
                id_records.append(newid)
            else:
                newid = ''.join([record.id, '_M'])
                id_records.append(newid)
            minhashes.append(E)
      
    simil = dict()
    for i, e in enumerate(minhashes):
        jac = list()
        for j, e2 in enumerate(minhashes):
            x = e.jaccard(minhashes[j])
            jac.append(x)
        simil[id_records[i]] = jac

    array = {k:np.array(v) for k, v in simil.items()}
    X = pd.DataFrame.from_dict(array, orient = 'index')
    sour_dist = pd.DataFrame.from_dict(simil)
    sour_path = '{outdir}/kmer'.format(outdir = outdir)
    sour_dist.to_csv(os.path.join(sour_path, '{prefix}_sourmash_distances.tsv'.format(prefix = prefix)), sep = '\t', index = False)
    
    plt.figure(figsize = (15, 10))
    sns.set(style = 'white', font_scale = 1.2, palette = 'Spectral')
    ax = sns.clustermap(X)
    sns.despine()
    ax.savefig(os.path.join(sour_path,'sourmash_clustermap.jpg'))
    plt.clf()
    







    
