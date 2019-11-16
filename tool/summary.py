#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 07:50:22 2019

@author: alerojo
"""

from Bio import SeqUtils, SeqIO
from Bio.Seq import Seq
import pandas as pd
import seaborn as sns
import os

# Function to extract the coordinates from the backbone file

def regions(reference, prefix, out):
    
    coordinates = '{out}/alignment/{genome_id}.backbone'.format(genome_id = prefix, out = out)
    
    # Parse backbone file 
    coordinates = pd.read_table(coordinates, sep = '\t')

    # Extract mapped regions coordinates
    mappedlocations = coordinates[(coordinates.seq1_leftend > 0) & (coordinates.seq1_rightend > 0)]
    mappedlocations = mappedlocations[['seq0_leftend','seq0_rightend']]
    mappedlocations = mappedlocations[(mappedlocations.seq0_leftend > 0) & (mappedlocations.seq0_rightend > 0)]

    # Extract unmapped regions coordinates
    unmappedlocations = coordinates[(coordinates.seq1_leftend == 0) & (coordinates.seq1_rightend == 0)]
    unmappedlocations = unmappedlocations[['seq0_leftend', 'seq0_rightend']]

    # Extract conflict regions coordinates
    conflictlocations = coordinates[(coordinates.seq0_leftend == 0) & (coordinates.seq0_rightend == 0)]
    conflictlocations = conflictlocations[['seq1_leftend', 'seq1_rightend']]
    
    return mappedlocations, unmappedlocations, conflictlocations
    
# Function to extract the regions from the reference and store them in dictionaries

def refextract(reference, mappedlocations, unmappedlocations, conflictlocations, prefix, flanking):
    
    # Parse reference FASTA file
    read = SeqIO.read(reference, format = 'fasta')
    
    # Create reference summary dictionary: GC content, length, number of mapped and unmapped regions
    #refstats_dict = dict()
    #refstats_dict = [{'GCContent': SeqUtils.GC(read.seq),
    #             'Length': len(str(read.seq)),
    #             'Unmapped': unmappedlocations.shape[0],
    #             'Mapped': mappedlocations.shape[0]}]
    
    # Extract mapped regions and store in a dictionary
    mappeddict = dict()
    idmap = list()
    
    for i in range(0, mappedlocations.shape[0]):
        start = mappedlocations.iloc[i,0]
        end = mappedlocations.iloc[i,1]
        mappeddict[i] = str(read.seq[start:end])
    
    for i in range(0, mappedlocations.shape[0]):
        start = mappedlocations.iloc[i,0]
        end = mappedlocations.iloc[i,1]
        header = (str(prefix),'_', str(start), ':', str(end))
        idmap.append(''.join(header))
    
    for i in range(0, len(mappeddict)):
        mappeddict[idmap[i]] = mappeddict.pop(i)
     
    # Extract unmapped regions and store in a dictionary
    unmappeddict = dict()
    idunmap = list()
    
    for i in range(0, unmappedlocations.shape[0]):
        start = unmappedlocations.iloc[i,0]
        end = unmappedlocations.iloc[i,1]
        unmappeddict[i] = str(read.seq[start-flanking:end+flanking])
    
    for i in range(0, unmappedlocations.shape[0]):
        start = unmappedlocations.iloc[i,0]
        end = unmappedlocations.iloc[i,1]
        header = (str(prefix),'_', str(start-flanking), ':', str(end+flanking))
        idunmap.append(''.join(header))

    for i in range(0, len(unmappeddict)):
        unmappeddict[idunmap[i]] = unmappeddict.pop(i)
    
    # Extract conflict regions and store in a dictionary
    conflictdict = dict()
    idconflict = list()
    
    for i in range(0, conflictlocations.shape[0]):
       start = conflictlocations.iloc[i,0]
       end = conflictlocations.iloc[i,1]
       conflictdict[i] = str(read.seq[start:end])
    
    for i in range(0, conflictlocations.shape[0]):
        start = conflictlocations.iloc[i,0]
        end = conflictlocations.iloc[i,1]
        header = (str(prefix),'_', str(start), ':', str(end))
        idconflict.append(''.join(header))

    for i in range(0, len(conflictdict)):
        conflictdict[idconflict[i]] = conflictdict.pop(i)
    
    return mappeddict, unmappeddict, idunmap, conflictdict

def unmapsum(unmappeddict, idunmap):
    
    # Create GC content, length and amino acid residues list to store values for each unmapped region
    gc_unmap = list()
    len_unmap = list()
    amino = pd.DataFrame(columns = ['A', 'D','E', 'G','F', 'L', 'Y', 'C', 'W', 'P', 'H', 'Q','I', 'M', 'T', 'N', 'S', 'K', 'R', 'V'])
    
    # Calculate values for each unmapped sequence
    for seq in unmappeddict.values():
        gc_unmap.append(SeqUtils.GC(str(seq)))
        len_unmap.append(len(seq))
        dna = Seq(seq)
        
        # Count number of residues for all six frames in the unmapped region sequence
        A = 0
        D = 0
        E = 0
        G = 0
        F = 0
        L = 0
        Y = 0
        C = 0
        W = 0
        P = 0
        H = 0
        Q = 0
        I = 0
        M = 0
        T = 0
        N = 0
        S = 0
        K = 0
        R = 0
        V = 0
        
        dna_seqs = [dna, dna.reverse_complement()]
        for s in dna_seqs:
            for i in range(3):
                
                aa_seq = s[i:].translate(table = 11)
            
                A = (str(aa_seq).count('A'))/len(str(aa_seq)) + A
                D = (str(aa_seq).count('D')/len(str(aa_seq))) + D
                E = (str(aa_seq).count('E'))/len(str(aa_seq)) + E
                G = (str(aa_seq).count('G'))/len(str(aa_seq)) + G
                F = (str(aa_seq).count('F'))/len(str(aa_seq)) + F
                L = (str(aa_seq).count('L'))/len(str(aa_seq)) + L
                Y = (str(aa_seq).count('Y'))/len(str(aa_seq)) + Y
                C = (str(aa_seq).count('C'))/len(str(aa_seq)) + C
                W = (str(aa_seq).count('W'))/len(str(aa_seq)) + W
                P = (str(aa_seq).count('P'))/len(str(aa_seq)) + P
                H = (str(aa_seq).count('H'))/len(str(aa_seq)) + H
                Q = (str(aa_seq).count('Q'))/len(str(aa_seq)) + Q
                I = (str(aa_seq).count('I'))/len(str(aa_seq)) + I
                M = (str(aa_seq).count('M'))/len(str(aa_seq)) + M
                T = (str(aa_seq).count('T'))/len(str(aa_seq)) + T
                N = (str(aa_seq).count('N'))/len(str(aa_seq)) + N
                S = (str(aa_seq).count('S'))/len(str(aa_seq)) + S
                K = (str(aa_seq).count('K'))/len(str(aa_seq)) + K
                R = (str(aa_seq).count('R'))/len(str(aa_seq)) + R
                V = (str(aa_seq).count('V'))/len(str(aa_seq)) + V
        
        amino = amino.append({'A':A*100,
                          'D':D*100,
                          'E':E*100,
                          'G':G*100,
                          'F':F*100,
                          'L':L*100,
                          'Y':Y*100,
                          'C':C*100,
                          'W':W*100,
                          'P':P*100,
                          'H':H*100,
                          'Q':Q*100,
                          'I':I*100,
                          'M':M*100,
                          'T':T*100,
                          'N':N*100,
                          'S':S*100,
                          'K':K*100,
                          'R':R*100,
                          'V':V*100}, ignore_index = True)
    
    # Create unmapped region summary dataframe: Region, GC content, length and total amino acid frequency for all six reading frames 
    unmap_stats = pd.DataFrame(list(zip(idunmap, gc_unmap, len_unmap)), columns = ['Region', 'GCContent', 'Length'])
    unmap_stats = pd.concat([unmap_stats, amino], axis = 1)
    unmap_stats.reset_index(drop = True, inplace = True)
    unmap_stats.sort_index(inplace = True)
    
    return unmap_stats

def refstats(reference, mappeddict, unmappeddict):
    
    # Calculate fraction of genome that is (un)mapped
    length_map = 0
    for key, values in mappeddict.items():
        length_map = length_map + len(values)
    
    length_un = 0
    for key, values in unmappeddict.items():
        length_un = length_un + len(values)
    
    read = SeqIO.read(reference, format = 'fasta')
    refstats_dict = dict()
    refstats_dict = [{'GCContent': SeqUtils.GC(read.seq),
                     'Length': len(str(read.seq)),
                     '#MappedRegions': len(mappeddict),
                     '#UnmappedRegions': len(unmappeddict),
                     'FractionMapped': (length_map/len(str(read.seq)))*100,
                     'FractionUnmapped': (length_un/len(str(read.seq)))*100}]
    
    # Create reference summary dataframe
    refstats_t = pd.DataFrame.from_dict(refstats_dict)
    refstats_t.reset_index(drop = True, inplace = True)
    refstats_t.sort_index(inplace = True)
    
    # Create reference summary dictionary: GC content, length, number of mapped and unmapped regions
    #refstats_dict = dict()
    #refstats_dict = [{'GCContent': SeqUtils.GC(read.seq),
    #             'Length': len(str(read.seq)),
    #             'Unmapped': unmappedlocations.shape[0],
    #             'Mapped': mappedlocations.shape[0]}]
    return refstats_t
        
def output(mappeddict, unmappeddict, conflictdict, refstats, unmap_stats, prefix, out):
    
    # Write summary tables
    newpath = 'summary'
    os.makedirs(os.path.join(out,newpath))
    path_sum = '{out}/summary'.format(out = out)
    refstats.to_csv(os.path.join(path_sum,'{genome_id}_referencesummary.tsv'.format(genome_id = prefix)), sep = '\t', index = False)
    unmap_stats.to_csv(os.path.join(path_sum,'{genome_id}_unmapsummary.tsv'.format(genome_id = prefix)), sep = '\t', index = False)
    
    # Write mapped regions FASTA file
    with open(os.path.join(out,'{prefix}_mappedregions.fasta'.format(prefix = prefix)), 'w+') as fasta:
        for key, value in mappeddict.items():
            fasta.write('>' + key + '\n' + value + '\n')
    
    # Write unmapped regions FASTA file
    with open(os.path.join(out,'{prefix}_unmappedregions.fasta'.format(prefix = prefix)), 'w+') as fasta:
        for key, value in unmappeddict.items():
            fasta.write('>' + key + '\n' + value + '\n')
    
    # Write mapped regions FASTA file
    with open(os.path.join(out,'{prefix}_conflictregions.fasta'.format(prefix = prefix)), 'w+') as fasta:
        for key, value in conflictdict.items():
            fasta.write('>' + key + '\n' + value + '\n')

#def plot()


def extract_main(reference, prefix, flanking, out):
    
    mappedlocations, unmappedlocations, conflictlocations = regions(reference, prefix, out)
    mappeddict, unmappeddict, idunmap, conflictdict = refextract(reference, mappedlocations, unmappedlocations, conflictlocations, prefix, flanking)
    unmap_stats = unmapsum(unmappeddict, idunmap)
    refstats_t = refstats(reference, mappeddict, unmappeddict)
    output(mappeddict, unmappeddict, conflictdict, refstats_t, unmap_stats, prefix, out)



    
    
