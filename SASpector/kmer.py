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
    newpath = 'kmer'
    os.makedirs(os.path.join(outdir,newpath))
    
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







    
