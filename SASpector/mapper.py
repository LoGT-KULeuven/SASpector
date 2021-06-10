#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: alerojo, 0mician
"""
import logging
import os
import subprocess
import shlex
import shutil


""" mapper

This script allows to align the reference and the Illumina contigs sequences using
progressiveMauve algorithm from Mauve command. The first input file is the reference FASTA file and the
second input file is the Illumina contigs FASTA file.

"""

def union(reference, prefix, out):
    """ Wraps union function from EMBOSS to check if reference is multifasta & combining
    
    Parameters
    ----------
    reference : str
        The file location of the reference FASTA file
    prefix : str
        Name of the genome
    out : str
        Output directory  
    """
    fasta_count = len([1 for line in open(reference) if line.startswith(">")])
          
    if(fasta_count > 1):
        logging.info("Your reference contains %i contigs. We are concatenating them into %s before pursuing" % (fasta_count, "{prefix}_concatenated.fasta".format(prefix=prefix)))
        cmd = 'union -sequence {reference} -outseq {concatenated}'.format(
            reference = reference, concatenated = "{prefix}_concatenated.fasta".format(prefix=prefix))
        process = subprocess.Popen(shlex.split(cmd), stdout = subprocess.DEVNULL, stderr=subprocess.STDOUT)
        process.wait()
        return True
    else:
        return False

        
def mauve(reference, contigs, prefix, out):
    """ Wraps progressiveMauve command line and generates the alignment outputs with the backbone file in alignment subdirectory
    
    Parameters
    ----------
    reference : str
        The file location of the reference FASTA file
    contigs : str
        The file location of the Illumina contigs FASTA file
    prefix : str
        Name of the genome
    out : str
        Output directory
    
    """
    logging.info("Starting whole genome alignment")
    cmd = 'progressiveMauve {reference} {contigs} --output={prefix}.alignment --backbone-output={prefix}.backbone'.format(
    reference = reference, contigs = contigs, prefix = prefix, out = out)
    process = subprocess.Popen(shlex.split(cmd), stdout = subprocess.DEVNULL, stderr=subprocess.STDOUT)
    process.wait()
    
    newdir = 'alignment'
    os.makedirs(os.path.join(out,newdir))
    shutil.move('{prefix}.alignment'.format(prefix = prefix), '{out}/alignment/{prefix}.alignment'.format(prefix = prefix, out = out))
    shutil.move('{prefix}.alignment.bbcols'.format(prefix = prefix), '{out}/alignment/{prefix}.alignment.bbcols'.format(out = out, prefix = prefix))
    shutil.move('{prefix}.backbone'.format(prefix = prefix), '{out}/alignment/{prefix}.backbone'.format(out = out, prefix = prefix))
    shutil.move('{reference}.sslist'.format(reference = reference), '{out}/alignment/{reference}.sslist'.format(out = out, reference = reference))
    shutil.move('{contigs}.sslist'.format(contigs = contigs), '{out}/alignment/{contigs}.sslist'.format(out = out, contigs = contigs))

    logging.info("Whole genome alignment completed!")
    

