#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 19:02:15 2019

@author: 0mician
"""
import glob
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys

import pandas as pd

from Bio import SeqIO
from Bio import Entrez

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def mash_output_parsing(results_file):
    logging.info("Now parsing the results file from mash")
    df = pd.read_table(results_file, header=None)
    df.sort_values(0, ascending=False, inplace=True)
    accn = re.search('complete_genomes/(.+?)_genomic.fna.gz', df.iloc[0,4]).group(1)
    accn = accn.split("_", 2)[:2]
    accn = "_".join(accn)
    logging.info("Results parsed! Found genome {accn}".format(accn = accn))
    return accn

def download_genome(accn, newdir, outdir):
    logging.info("Downloading genome {accn} from NCBI".format(accn = accn))
    Entrez.email = "SASpector@omician.net"
    esearch_handle = Entrez.esearch(db="nuccore", term=accn)
    esearch_results = Entrez.read(esearch_handle)
    esearch_handle.close()
    nuccore_id = esearch_results['IdList'][0]
    efetch_handle = Entrez.efetch(db="nuccore", id=nuccore_id, rettype="fasta", retmode="text")
    out_handle = open(os.path.join(outdir, newdir, "{accn}.fasta".format(accn=accn)), "w")
    out_handle.write(efetch_handle.read())
    out_handle.close()
    logging.info("Genome downloaded and saved in {outdir}/mash/{accn}.fasta".format(outdir=outdir, accn=accn))
    return "{outdir}/mash/{accn}.fasta".format(outdir=outdir, accn=accn)

def calculate_anib(draft, reference_selected, prefix, outdir):
    logging.info("Let's assess the similarity of the reference selected using the ANI metrics")
    logging.info("SASpector will now run pyani using the blast method")

    try:
        os.makedirs(os.path.join(outdir, "ANI/ANIb_input"))
        shutil.copy(draft, '{outdir}/ANI/ANIb_input/'.format(outdir = outdir))
        shutil.copy(reference_selected, '{outdir}/ANI/ANIb_input/'.format(outdir = outdir))

        cmd = 'average_nucleotide_identity.py -m ANIb -i {outdir}/ANI/ANIb_input -o {outdir}/ANI/ANIb_results'.format(outdir = outdir)
        process = subprocess.run(shlex.split(cmd), stderr = subprocess.DEVNULL)
        logging.info("ANIb process completed, now parsing the results")

    except IOError:
        logging.error("There was an issue during the ANI analysis")
        logging.error("Exiting SASpector")
        sys.exit()
    
    identity_perc_file = os.path.join(outdir, "ANI/ANIb_results/ANIb_percentage_identity.tab")
    df = pd.read_table(identity_perc_file)
    avg_anib = ((df.iloc[1,1] + df.iloc[0, 2])/2.0)*100
    #df.sort_values(0, ascending=False, inplace=True)
    logging.info("ANI value: {avg_anib}".format(avg_anib = avg_anib))
    
def sgmash(draft, refseq_msh, prefix, outdir):
    """ Wraps mash and selects a closed assembly from refseq
    
    Parameters
    ----------
    draft : str
        The file location of the draft assembly FASTA file
    refseq_msh: str
        sketch file of the refseq complete genomes
    prefix: str
        Name of the genome
    outdir: str
        Output directory
    
    """
    logging.info("Running genome selection with mash")
    logging.info("This feature is experimental!")
    logging.info("It is possible that a genome is autoselected which does not represent your draft assembly.")

    newdir = 'mash'
    os.makedirs(os.path.join(outdir, newdir))
    cmd = 'mash screen {refseq_msh} {draft}'.format(refseq_msh = refseq_msh, draft = draft)
    path_results = os.path.join(outdir, newdir, "results.txt")
    file_results = open(path_results, "w")
    process = subprocess.run(shlex.split(cmd), stdout = file_results, stderr = subprocess.DEVNULL)
    logging.info("Mash screen process completed")

    if(is_non_zero_file(path_results)):
        accn = mash_output_parsing(path_results)
        reference_selected = download_genome(accn, newdir, outdir)
    else:
        logging.info("Something went wrong with the selection.")
        sys.exit("Program stopped!")

    calculate_anib(draft, reference_selected, prefix, outdir)
    
    return reference_selected
