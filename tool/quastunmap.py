#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:46:34 2019

@author: alerojo
"""

import subprocess
import shlex
import progressbar 
import time

def quast(reference, outdir, prefix):
    bar = progressbar.ProgressBar(widgets = ['Running QUAST: ', progressbar.Bar(), '(', progressbar.ETA(),')'])
    for i in bar(range(1)):
        cmd = 'quast.py {outdir}/{prefix}_unmappedregions.fasta -r {reference} -o {outdir}/quast'.format(
                outdir = outdir, reference = reference, prefix = prefix)
        process = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE)
        while process.poll() is None:
            l = process.stdout.readline()
    time.sleep(0.02)
    

    
