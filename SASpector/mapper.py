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
    newdir = 'alignment'
    os.makedirs(os.path.join(out,newdir))
    shutil.copy(reference, '{out}/alignment/'.format(out = out))
    shutil.copy(contigs, '{out}/alignment/'.format(out = out))

    reference_file = os.path.join(out, newdir, os.path.basename(reference))
    contigs_file = os.path.join(out, newdir, os.path.basename(contigs))
    alignment_file = os.path.join(out, newdir, '{prefix}.alignment'.format(prefix=prefix))
    backbone_file = os.path.join(out, newdir, '{prefix}.backbone'.format(prefix=prefix))

    cmd = 'progressiveMauve {reference_file} {contigs_file} --output={alignment_file} --backbone-output={backbone_file}'.format(
        reference_file = reference_file, contigs_file = contigs_file, alignment_file=alignment_file, backbone_file=backbone_file)
    process = subprocess.Popen(shlex.split(cmd), stdout = subprocess.DEVNULL, stderr=subprocess.STDOUT)
    process.wait()
    
    logging.info("Whole genome alignment completed!")
    

