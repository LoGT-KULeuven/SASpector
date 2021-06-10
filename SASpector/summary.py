#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 07:50:22 2019

@author: alerojo, 0mician
"""

from Bio import SeqUtils, SeqIO, BiopythonWarning
from Bio.Seq import Seq
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import logging
import os
import time
import warnings


""" summary

This script takes the output of progressiveMauve alignment (the backbone file) and extracts the mapped regions, conflict regions and missing regions from
the reference genome. It parses the backbone file and generates a FASTA file with the extracted regions. Additionally, it generates basic 
sequence statistics such as length, GC content and amino acid residues frequency percentages for both mapped and unmapped regions. 

"""

# Function to extract the coordinates from the backbone file

def regions(prefix, out):
    """ Parses the backbone file and extracts the coordinates of the mapped regions, missing regions and conflict regions.
        Conflict regions are regions that don't map perfectly to the reference due to gaps or indels.
    
    Parameters
    ----------
    prefix: str
        Name of genome
    out: str
        Output directory
    
    Returns
    -------
    mappedlocations: dataframe
        Coordinates of the mapped regions in the reference sequence
    unmappedlocations: dataframe
        Coordinates of the unmapped regions in the reference sequence
    conflictlocations: dataframe
        Coordinates of the conflict regions in the reference sequence
    reverselocations: dataframe
        Coordinates of the mapped reverse complementary regions in the reference sequence
        
    """
    coordinates = '{out}/alignment/{genome_id}.backbone'.format(genome_id = prefix, out = out)
    
    # Parse backbone file 
    coordinates = pd.read_csv(coordinates, sep = '\t')

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
    
    # Extract mapped regions coordinates
    mappedlocations = coordinates[(coordinates.seq1_leftend > 0) & (coordinates.seq1_rightend > 0)]
    mappedlocations = mappedlocations[['seq0_leftend','seq0_rightend']]
    mappedlocations = mappedlocations[(mappedlocations.seq0_leftend > 0) & (mappedlocations.seq0_rightend > 0)]
    
    # Extract reverse coordinates
    reverselocations = coordinates[(coordinates.seq1_leftend > 0) & (coordinates.seq1_rightend > 0)]
    reverselocations = reverselocations[['seq0_leftend','seq0_rightend']]
    reverselocations = reverselocations[(reverselocations.seq0_leftend < 0) & (reverselocations.seq0_rightend < 0)]
    
    return mappedlocations, unmappedlocations, conflictlocations, reverselocations
    
# Function to extract the regions (+ flanks) from the reference and store them in dictionaries

def refextract(reference, mappedlocations, unmappedlocations, conflictlocations, prefix, flanking):
    """ Extracts the regions from the reference FASTA file
    
    Parameters
    ----------
    reference: str
        The file location of the reference FASTA file
    prefix: str
        Name of the genome
    flanking: int
        Length of flanking regions [Default = 0 bp]
    mappedlocations: dataframe
        Coordinates of the mapped regions in the reference sequence
    unmappedlocations: dataframe
        Coordinates of the unmapped regions in the reference sequence
    conflictlocations: dataframe
        Coordinates of the conflict regions in the reference sequence
    reverselocations: dataframe
        Coordinates of the mapped reverse complement regions in the reference sequence
    
    Returns
    -------
    mappeddict: dict
        Dictionary of the coordinates and sequences of the mapped regions
    unmappeddict: dict
        Dictionary of the coordinates and sequences of the unmapped regions
    idunmap: list
        List of the coordinates of all unmapped regions as strings
    conflictdict: dict
        Dictionary of the coordinates and sequences of the conflict regions
    
    """
    # Parse reference FASTA file
    read = SeqIO.read(reference, format = 'fasta')
        
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
        header = (str(prefix), '_',str(start), ':', str(end))
        idmap.append(''.join(header))
    
    for i in range(0, len(mappeddict)):
        mappeddict[idmap[i]] = mappeddict.pop(i)
     
    # Extract unmapped regions and store in a dictionary
    unmappeddict = dict()
    idunmap = list()
    
    for i in range(0, unmappedlocations.shape[0]):
        start = unmappedlocations.iloc[i,0]
        end = unmappedlocations.iloc[i,1]
        if len(str(read.seq[start-flanking:end+flanking])) > 100:
            unmappeddict[i] = str(read.seq[start-flanking:end+flanking])
    unmappeddict = {i: v for i, v in enumerate(unmappeddict.values())}
    
    for i in range(0, unmappedlocations.shape[0]):
        start = unmappedlocations.iloc[i,0]
        end = unmappedlocations.iloc[i,1]
        if len(str(read.seq[start-flanking:end+flanking])) > 100:
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
        header = (str(prefix), '_', str(start), ':', str(end))
        idconflict.append(''.join(header))

    for i in range(0, len(conflictdict)):
        conflictdict[idconflict[i]] = conflictdict.pop(i)
    

    return mappeddict, unmappeddict, idunmap, conflictdict

def unmapsum(unmappeddict, idunmap):
    """ Generates summary statistics for the extracted unmapped regions
    
    Parameters
    ----------
    unmappeddict: dict
        Dictionary of the coordinates and sequences of the unmapped regions
    idunmap: list
        List of the coordinates for all unmapped regions as strings
    
    Returns
    -------
    unmap_stats: dataframe
        Table containing the coordinates, GC content, length and amino acid residues frequency percentages for each missing region
    
    """
    
    # Create GC content, length and amino acid residues list to store values for each unmapped region
    gc_unmap = list()
    len_unmap = list()
    amino = pd.DataFrame(columns = ['A', 'D','E', 'G','F', 'L', 'Y', 'C', 'W', 'P', 'H', 'Q','I', 'M', 'T', 'N', 'S', 'K', 'R', 'V'])
    
    # Calculate values for each unmapped sequence
    for seq in unmappeddict.values():
        gc_unmap.append(SeqUtils.GC(str(seq)))
        len_unmap.append(len(seq))
        dna = Seq(seq)
        dna_seq = [dna, dna.reverse_complement()]
        codes = list()
        for s in dna_seq:
            for frame in range(3):
                pro = s[frame:].translate(table = 11)
                codes.append(pro._data)
                
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
        Stop = 0
        len_seq = 0
            
        for pro_seq in codes:
                
            # Counts
            A = (str(pro_seq).count('A')) + A
            D = (str(pro_seq).count('D')) + D
            E = (str(pro_seq).count('E')) + E
            G = (str(pro_seq).count('G')) + G
            F = (str(pro_seq).count('F')) + F
            L = (str(pro_seq).count('L')) + L
            Y = (str(pro_seq).count('Y')) + Y
            C = (str(pro_seq).count('C')) + C
            W = (str(pro_seq).count('W')) + W
            P = (str(pro_seq).count('P')) + P
            H = (str(pro_seq).count('H')) + H
            Q = (str(pro_seq).count('Q')) + Q
            I = (str(pro_seq).count('I')) + I
            M = (str(pro_seq).count('M')) + M
            T = (str(pro_seq).count('T')) + T
            N = (str(pro_seq).count('N')) + N
            S = (str(pro_seq).count('S')) + S
            K = (str(pro_seq).count('K')) + K
            R = (str(pro_seq).count('R')) + R
            V = (str(pro_seq).count('V')) + V
            Stop = (str(pro_seq).count('*')) + Stop
            len_seq = len(pro_seq) + len_seq
        
        A = A/len_seq
        D = D/len_seq
        E = E/len_seq
        G = G/len_seq
        F = F/len_seq
        L = L/len_seq
        Y = Y/len_seq
        C = C/len_seq
        W = W/len_seq
        P = P/len_seq
        H = H/len_seq
        Q = Q/len_seq
        I = I/len_seq
        M = M/len_seq
        T = T/len_seq
        N = N/len_seq
        S = S/len_seq
        K = K/len_seq
        R = R/len_seq
        V = V/len_seq
        Stop = Stop/len_seq
        
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
                          'V':V*100,
                          'Stop':Stop*100}, ignore_index = True)
    codes.clear()
                
    
    # Create unmapped region summary dataframe: Region, GC content, length and total amino acid frequency for all six reading frames 
    unmap_stats = pd.DataFrame(list(zip(idunmap, gc_unmap, len_unmap)), columns = ['Region', 'GCContent', 'Length'])
    unmap_stats = pd.concat([unmap_stats, amino], axis = 1)
    unmap_stats.reset_index(drop = True, inplace = True)
    unmap_stats.sort_index(inplace = True)
    
    return unmap_stats

def refstats(reference, mappedlocations, unmappedlocations, conflictlocations, reverselocations, unmappeddict):
    """Generates summary statistics for reference genome based on the mapped, unmapped and conflict regions
    
    Parameters
    ----------
    reference: str
        The file location of the reference FASTA file
    mappedlocations: dataframe
        Coordinates of the mapped regions in the reference sequence
    unmappedlocations: dataframe
         Coordinates of the unmapped regions in the reference sequence
    conflictlocations: dataframe
        Coordinates of the conflict regions in the reference sequence
    reverselocations: dataframe
        Coordinates of the mapped reverse complement regions in the reference sequence
    unmappeddict: dataframe
        Dictionary of the coordinates and sequences of the unmapped regions
    
    Returns
    -------
    refstats_t: dataframe
        Table containing the reference summary statistics
    
        
    """
    # Calculate genome fraction
    sum_map = 0
    for i in range(0, mappedlocations.shape[0]):
        sum_map = sum_map + abs(mappedlocations.iloc[i,1] - mappedlocations.iloc[i,0])
    
    sum_confl = 0
    for i in range(0, conflictlocations.shape[0]):
        sum_confl = sum_confl + abs(conflictlocations.iloc[i,1] - conflictlocations.iloc[i,0])
        
    sum_rev = 0
    for i in range(0, reverselocations.shape[0]):
        sum_rev = sum_rev + abs(reverselocations.iloc[i,1] - reverselocations.iloc[i,0])
        
    total_map = sum_map + sum_confl + sum_rev
    
    sum_unmap = 0
    for i in range(0, unmappedlocations.shape[0]):
        sum_unmap = sum_unmap + abs(unmappedlocations.iloc[i,1] - unmappedlocations.iloc[i,0])
    
    read = SeqIO.read(reference, format = 'fasta')
    refstats_dict = dict()
    refstats_dict = [{'GCContent': SeqUtils.GC(read.seq),
                     'Length': len(str(read.seq)),
                     'NumberMappedRegions': mappedlocations.shape[0] + reverselocations.shape[0] + conflictlocations.shape[0],
                     'NumberUnmappedRegions': unmappedlocations.shape[0],
                     'FilteredUnmappedRegions': len(unmappeddict),
                     'FractionMapped': (total_map/len(str(read.seq)))*100,
                     'FractionUnmapped': (sum_unmap/len(str(read.seq)))*100}]
    
    # Create reference summary dataframe
    refstats_t = pd.DataFrame.from_dict(refstats_dict)
    refstats_t.reset_index(drop = True, inplace = True)
    refstats_t.sort_index(inplace = True)
    
    return refstats_t
        
def output(mappeddict, unmappeddict, conflictdict, refstats, unmap_stats, prefix, out):
    """Generates the FASTA files and summary tables (csv files) for the mapped regions, missing regions and conflict regions
    
    Parameters
    ----------
    mappeddict: dict
        Dictionary of the coordinates and sequences of the mapped regions
    conflictddict: dict
        Dictionary of the coordinates and sequences of the conflict regions
    unmappeddict: dict
        Dictionary of the coordinates and sequences of the unmapped regions
    refstats: dataframe
        Table containing the summary statistics of the reference
    unmap_stats: dataframe
        Table containing the summary statistics of each missing region
    prefix: str
        Name of genome
    out: str
        Output directory
        
    """
    
    path_sum = '{out}'.format(out = out)
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

def plot(unmappeddict, unmap_stats, out):
    """ Generates boxplot, distribution plots and join plots from the missing regions summary statistics.
        They are saved in the output directory as jpg images.
    
    Parameters
    ----------
    unmappeddict: dict
        Dictionary of the coordinates and sequences of the unmapped regions
    unmap_stats: dataframe
        Table containing the unmapped regions summary statistics
    out: str
        Output directory
        
    """
    
    gc_content = list()
    regions_length = list()
    for key, values in unmappeddict.items():
        gc_content.append(SeqUtils.GC(values))
        regions_length.append(len(values))
    
    plt.figure(figsize = (10,10))
    sns.set(style = 'white', font_scale = 2)
    fig_joint = sns.jointplot(x = regions_length, y = gc_content, kind = 'hex', height = 7)
    fig_joint.set_axis_labels(xlabel = 'Length', ylabel = 'GC Content')
    fig_joint.savefig(os.path.join(out, 'gc_length_joint_missing.jpg'))
    plt.clf()

    plt.figure(figsize = (15,10))
    sns.set(style = 'white', font_scale = 1.3)
    fig_gc = sns.displot(gc_content, kind = "hist", rug = False, color = 'red')
    fig_gc.set(xlabel = 'GC Content')
    fig_gc.fig.suptitle('Distribution of GC Content')
    sns.despine()
    fig_gc.savefig(os.path.join(out, 'gc_content_missing.jpg'))
    plt.clf()

    plt.figure(figsize = (15,10))
    sns.set(style = 'white', font_scale = 1.3)
    fig_length = sns.displot(regions_length, kind = "hist", rug = False, color = 'green')
    fig_length.set(xlabel = 'Length')
    fig_length.fig.suptitle('Distribution of Length')
    sns.despine()
    fig_length.savefig(os.path.join(out, 'length_missing.jpg'))
    plt.clf()
    
    plt.figure(figsize = (15, 10))
    sns.set(style = 'white', font_scale = 1.2)
    codons = sns.boxplot(data = unmap_stats.iloc[:, 3:24], palette = 'Spectral')
    codons.set(xlabel = 'Translated Codons', ylabel = 'Mean Percentage per Frame (%)')
    sns.despine()
    codons.figure.savefig(os.path.join(out, 'codons_missing.jpg'))
    plt.clf()


def extract_main(reference, prefix, flanking, out):
    """ Main function of this script
    
    Parameters
    ----------
    reference: str
        The file location of the reference FASTA file
    prefix: str
        Name of genome
    flanking: int
        Length of the flanking regions [Default = 0 bp]
    out: str
        Output directory
    
    Returns
    -------
    mappedlocations: dataframe
        Coordinates of the mapped regions
    unmappedlocations: dataframe
        Coordinates of the unmapped regions
    conflictlocations: dataframe
        Coordinates of the conflict regions
    reverselocations: dataframe
        Coordinates of the mapped reverse complementary regions
        
    """
    logging.info("Analysis of the whole genome alignment and extraction of regions of interest")
    warnings.simplefilter('ignore', BiopythonWarning)
    mappedlocations, unmappedlocations, conflictlocations, reverselocations = regions(prefix, out)
    mappeddict, unmappeddict, idunmap, conflictdict = refextract(reference, mappedlocations, unmappedlocations, conflictlocations, prefix, flanking)
    unmap_stats = unmapsum(unmappeddict, idunmap)
    refstats_t = refstats(reference, mappedlocations, unmappedlocations, conflictlocations, reverselocations, unmappeddict)
    plot(unmappeddict, unmap_stats, out)
    time.sleep(0.02)
    output(mappeddict, unmappeddict, conflictdict, refstats_t, unmap_stats, prefix, out)
    return mappedlocations, unmappedlocations, conflictlocations, reverselocations


