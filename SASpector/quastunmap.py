#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:46:34 2019

@author: alerojo, 0mician
"""
import logging
import subprocess
import shlex

""" quastunmap

This script allows to run QUAST to perform a genome quality assessment of the missing regions with the reference Hybrid assembly.
The input files are the missing regions FASTA file and the reference FASTA file.

"""

def quast(reference, outdir, prefix):
    """ Wraps QUAST and generates a genome assessment report, including a genome viewer
    
    Parameters
    ----------
    reference: str
        The file location of the reference FASTA file
    prefix: str
        Name of the genome
    outdir: str
        Output directory
    
    """
    logging.info("Running quast of draft genome against the reference")
    cmd = 'quast.py {outdir}/{prefix}_unmappedregions.fasta -r {reference} -o {outdir}/quast'.format(
        outdir = outdir, reference = reference, prefix = prefix)
    process = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE)
    process.wait()
    logging.info("QUAST completed")
    

    
