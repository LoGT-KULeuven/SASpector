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

import argparse
import os



def main():
    parser = argparse.ArgumentParser(prog = 'SASpector - Short-read Assembly inSpector', description = '')
    parser.add_argument('reference', type = str, help = 'Hybrid assembly FASTA file as reference genome')
    parser.add_argument('contigs', type = str, help = 'Illumina FASTA file as contigs/draft genome')
    parser.add_argument('-p', '--prefix', type = str, help = 'Genome ID')
    parser.add_argument('-dir', '--outdir', help = 'Output directory')
    parser.add_argument('-f', '--flanking', nargs = '?' ,const = 'flanking', type = int, help = 'Add flanking regions [Default = 0]', default = 0)
    parser.add_argument('-k', '--kmers', nargs = '?' ,const = 'kmers', type = int, help = 'Calculate kmer frequencies', default = 0)
    parser.add_argument('-q','--quast', help = 'Run QUAST for unmapped regions against reference assembly', action = 'store_true')

    args = parser.parse_args()
    
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    
    mauve(args.reference, args.contigs, args.prefix, args.outdir)
    extract_main(args.reference, args.prefix, args.flanking, args.outdir)
    prokka(args.prefix, args.outdir)
    blast(args.outdir, args.prefix)
    if args.kmers:
        kmer(args.kmers, args.prefix, args.outdir)
    if args.quast is True:
         quast(args.reference, args.outdir, args.prefix)
    
    print('Done!')

if __name__ == '__main__':
    main()

