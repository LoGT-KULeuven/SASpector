#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 11:46:33 2019

@author: alerojo, 0mician
"""       
import logging
import os
import subprocess
import shlex
import sys

from Bio.Blast.Applications import NcbiblastxCommandline as blastx

""" gene_predict

This script predicts genes that are in the missing regions of the reference FASTA file using Prokka and performs
a blastx to search hits across several bacterial species from an in-house protein database using BLAST+. The input file for Prokka
is the missing regions FASTA file generated in the summary script. The input files for blastx are the predicted genes DNA sequences FASTA file
and the bacterial protein sequences in FASTA file.

"""

def prokka(prefix, outdir):
    """ Wraps Prokka and generates the predicted genes per missing region with their annotation
    
    Parameters
    ----------
    prefix: str
        Name of the genome
    outdir: str
        Output directory
        
    """
    logging.info("Starting annotation with prokka")
    
    try:
        cmd = 'prokka --outdir {outdir}/genesprediction '\
              '--prefix {prefix}.predictedgenes '\
              '{outdir}/{prefix}_unmappedregions.fasta'.format(prefix = prefix, outdir = outdir)
        process = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        output,error = process.communicate()
        if process.returncode != 0:
            logging.error("Something went wrong with the prokka annotation")
            logging.error("Prokka exited with code {returncode}".format(returncode = process.returncode))
            logging.error("And the following error message {error}".format(error = error.strip()))
            sys.exit()
    except OSError as e:
        logging.error("OSError > {errno}".format(errno = e.errno))
        logging.error("OSError > {strerror}".format(strerror = e.strerror))
        logging.error("OSError > {filename}".format(filename = e.filename))
        sys.exit()
    except:
        logging.error("Exiting SASpector {exc_info}".format(exc_info = sys.exc_info()[0]))
        sys.exit()
        
    logging.info("Annotation process completed")
            
def blast(outdir, prefix, proteindb):
    """ Wraps BLAST+ and performs blastx with output Prokka nucleotide FASTA file and SASpector protein FASTA file
    
    Parameters
    ----------
    prefix: str
        Name of the genome
    outdir: str
        Output directory
    proteindb: str
        The file location of the protein database FASTA file
    
    """
    logging.info("Blasting of genes found in missing regions against "\
                 "protein fasta file provided (or default file if none provided")
    try:
        cline = blastx(cmd = 'blastx', out = '{outdir}/{prefix}_blastxresults.tsv'.format(prefix = prefix, outdir = outdir), evalue = 0.001,  outfmt = '6 qacc sacc stitle qlen slen qstart qend sstart send sstrand length nident mismatch positive evalue, stderr = cline()stdout, stderr = cline() send pident evalue qcovs', query = '{outdir}/genesprediction/{prefix}.predictedgenes.fna'.format(outdir = outdir, prefix = prefix), subject = proteindb, max_target_seqs = 5)
        stdout, stderr = cline()
    except:
        logging.error("BLAST error, maybe disable the --proteindb option?)")
        logging.error("Exiting SASpector")
        sys.exit()
        
    logging.info("Blast completed")

    
