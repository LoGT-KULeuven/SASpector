#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 19:02:15 2019

@author: alerojo, 0mician
"""

from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import os
import logging
import shlex
import shutil
import subprocess
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
    logging.info("Running k-mer analysis")    
    newpath = 'kmer'
    os.makedirs(os.path.join(outdir,newpath))
    
    # Define the unmapped regions FASTA file
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

        if all(value == 1 for value in kmers.values()): # no need to save/create figure in that case
            continue
        
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
    logging.info("K-mer analysis completed")
        
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
    logging.info("Running clustermap analysis with Sourmash")
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
    
    import warnings # locally deactivating warning about the use of fastcluster (not applicable)
    warnings.simplefilter("ignore")
    ax = sns.clustermap(X)

    sns.despine()
    ax.savefig(os.path.join(sour_path,'sourmash_clustermap.jpg'))
    plt.clf()
    logging.info("Clustermap analysis complete")







    
