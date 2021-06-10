#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 19:02:15 2019

@author: alerojo, 0mician
"""
import glob
import logging
import os
import shlex
import shutil
import subprocess



def trf(prefix,outdir):
    """ Wraps TRF and generates tandem repeats reports for the missing regions
    
    Parameters
    ----------
    prefix: str
        Name of the genome
    outdir: str
        Output directory
    
    """
    logging.info("Running the tandem repeats detection with trf")
    cmd = 'trf {outdir}/{prefix}_unmappedregions.fasta 2 5 7 80 10 50 2000'.format(prefix = prefix, outdir = outdir)
    process = subprocess.run(shlex.split(cmd), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    newdir = 'trf'
    os.makedirs(os.path.join(outdir, newdir))
    dest = os.path.join(outdir,newdir)
    for html in glob.glob('*.html'):
        shutil.move(html,dest)
    logging.info("Tandem repeat annotation completed")
