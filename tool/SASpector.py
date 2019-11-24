#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 06:02:38 2019

@author: alerojo
"""

from mapper import mauve
from summary import extract_main
from gene_predict import prokka, blast
from kmer import kmer
from quastunmap import quast
from coverage import cvg_main

import argparse
import os



def main():
    parser = argparse.ArgumentParser(prog = 'SASpector - Short-read Assembly inSpector', description = 'Exctract and analyse missing regions from short-read assemblies.')
    parser.add_argument('reference', type = str, help = 'Reference genome FASTA file, e.g. from hybrid assembly.\n If the file contains multiple seqences, only the first one is used, so make sure to concatenate if needed.')
    parser.add_argument('contigs', type = str, help = 'Contigs FASTA file from short-read assembly')
    parser.add_argument('-p', '--prefix', type = str, help = 'Genome ID')
    parser.add_argument('-dir', '--outdir', help = 'Output directory')
    parser.add_argument('-f', '--flanking', nargs = '?', const = 'flanking', type = int, help = 'Add flanking regions [Default = 0 bp]', default = 0)
    parser.add_argument('-k', '--kmers', nargs = '?', const = 'kmers', type = int, help = 'Calculate kmer frequencies (provide k)', default = 0)
    parser.add_argument('-q','--quast', help = 'Run QUAST for unmapped regions against reference assembly', action = 'store_true')
    parser.add_argument('-c', '--coverage', nargs='?', const='coverage', metavar='BAMFILE', type = str, help = 'Run SAMtools bedcov to look at short-read coverage in the missing regions. Needs alignment of reads to the reference genome in bam format')

    args = parser.parse_args()
    
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    
    mauve(args.reference, args.contigs, args.prefix, args.outdir)
    mappedlocations, unmappedlocations, conflictlocations, reverselocations = extract_main(args.reference, args.prefix, args.flanking, args.outdir)
    prokka(args.prefix, args.outdir)
    blast(args.outdir, args.prefix)
    if args.kmers:
        kmer(args.kmers, args.prefix, args.outdir)
    if args.quast is True:
        quast(args.reference, args.outdir, args.prefix)
    if args.coverage:
        cvg_main(mappedlocations, conflictlocations, args.coverage, args.reference, args.outdir, args.prefix)
    
    print('Done!')

if __name__ == '__main__':
    main()

