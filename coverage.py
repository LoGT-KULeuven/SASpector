# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 21:58:29 2019

@author: emmav
"""

import pandas as pd
import csv
from Bio import SeqIO
import subprocess
import os
import progressbar
import seaborn as sns
import matplotlib.pyplot as plt

""" coverage

This script calculate the average coverage of the mapped and unmapped regions. Using SAMtools, it calculates
the coverage for both regions and generates a summary table including the coordinates for each region, total base-pair length
and average depth. Additionally, it generates a barplot of the coverage for each mapped and unmapped region. The input file is
a BAM file from the reference and the Illumina reads.

"""

# bwa mem {inputRef} {input.R1} {input.R2} | samtools view -b -o -- | bwa index-a bwtsw {inputREF}

def make_bed(mappedlocations, conflictlocations, reference, outdir, prefix):
    """ Generates a bed file from the mapped locations, unmapped locations and conflic locations
    
    Parameters
    ----------
    mappedlocations: dataframe
        Coordinates of the mapped regions
    conflictlocations: dataframe
       Coordinates of the conflict regions
    reference: str
        The file location of the reference FASTA file
    outdir: str
        Output directory
    prefix: str
        Name of the genome
        
    """
    for seq in SeqIO.parse(reference, "fasta"):
        ID = seq.id.split(' ')[0]

    #Only the filtered unmapped regions from summary file
    unmap = []
    unmapsum = '{outdir}/{prefix}_unmapsummary.tsv'.format(outdir = outdir, prefix = prefix)
    regions = pd.read_csv(unmapsum, sep = '\t')
    regions = regions['Region'].values.tolist()
    for region in regions:
        start = region.split('_')[1].split(':')[0]
        end = region.split('_')[1].split(':')[1]
        dictunmap = {'id': ID, 'start': start, 'end': end}
        unmap.append(dictunmap)

    #Mapped regions: from df generated in summary script
    mapped = []
    mappedlocations.columns = ['start', 'end']
    conflictlocations.columns = ['start', 'end']
    mappedloc = mappedlocations.append(conflictlocations)
    mappedloc = mappedloc.sort_values(by=['start'])
    for i in mappedloc.index.values.tolist():
        dictmap = {'id': ID, 'start': mappedloc.loc[i, 'start'], 'end': mappedloc.loc[i, 'end']}
        mapped.append(dictmap)

    bedunmap = '{outdir}/coverage/{prefix}_unmappedregions.bed'.format(outdir = outdir, prefix = prefix)
    bedmap = '{outdir}/coverage/{prefix}_mappedregions.bed'.format(outdir = outdir, prefix = prefix)

    with open(bedunmap, "w") as outunmap, open(bedmap, "w") as outmap:
        writermap = csv.DictWriter(outmap, fieldnames=['id', 'start', 'end'], delimiter='\t')
        for region in mapped:
            writermap.writerow(region)
        writerunmap = csv.DictWriter(outunmap, fieldnames=['id', 'start', 'end'], delimiter='\t')
        for region in unmap:
            writerunmap.writerow(region)


def sam(bamfile, outdir, prefix): 
    """Wraps SAMtools to generate a sorted BAM file and the coverage for both mapped and unmapped regions
    
    Parameters
    ----------
    bamfile: str
        The file location of the BAM file
    outdir: str
        Output directory
    prefix: str
        Name of genome
        
    """
    sort = 'samtools sort -o {outdir}/coverage/{prefix}.sorted.bam {bam}'.format(outdir = outdir, prefix = prefix, bam = bamfile)
    index = 'samtools index {outdir}/coverage/{prefix}.sorted.bam'.format(outdir = outdir, prefix = prefix)
    bedcovu = 'samtools bedcov {outdir}/coverage/{prefix}_unmappedregions.bed {outdir}/coverage/{prefix}.sorted.bam > {outdir}/coverage/{prefix}_unmapcvg.tsv'.format(outdir = outdir, prefix = prefix)
    bedcovm = 'samtools bedcov {outdir}/coverage/{prefix}_mappedregions.bed {outdir}/coverage/{prefix}.sorted.bam > {outdir}/coverage/{prefix}_mapcvg.tsv'.format(outdir = outdir, prefix = prefix)

    pipe1 = subprocess.Popen(sort, shell = True, stdout = subprocess.PIPE)
    pipe1.wait()
    pipe2 = subprocess.Popen(index, shell = True, stdin = pipe1.stdout, stdout = subprocess.PIPE)
    pipe2.wait()
    pipe3 = subprocess.Popen(bedcovu, shell = True, stdin = pipe2.stdout, stdout = subprocess.PIPE)
    pipe3.wait()
    pipe4 = subprocess.Popen(bedcovm, shell = True, stdin = pipe3.stdout, stdout = subprocess.PIPE)
    pipe4.wait()
    
def output(outdir, prefix): 
    """Writes the coverage statistics for each mapped and unmapped regions, and generates a barplot of the coverage for each region
    
    Parameters
    ----------
    outdir: str
        Output directory
    prefix: str
        Name of genome
    """
    bedcovu = '{outdir}/coverage/{prefix}_unmapcvg.tsv'.format(prefix = prefix, outdir = outdir)
    bedcovm = '{outdir}/coverage/{prefix}_mapcvg.tsv'.format(prefix = prefix, outdir = outdir)

    cvgu = pd.read_csv(bedcovu, header=None, sep = '\t')
    cvgu['perbaseavg'] = cvgu[3]/(cvgu[2]-cvgu[1])
    
    #Compare with mapped regions
    cvgm = pd.read_csv(bedcovm, header=None, sep = '\t')
    cvgm['perbaseavg'] = cvgm[3]/(cvgm[2]-cvgm[1])

    cvgu['flag']='unmapped'
    cvgm['flag']='mapped'
    cvg=cvgu.append(cvgm, ignore_index = True)
    cvg.columns=['seq','start', 'end', 'total_perbase_depth', 'avg_perbase_depth', 'flag']
   
    id_region = list()
    for i in range(0,cvg.shape[0]):
        start = str(cvg.iloc[i,1])
        end = str(cvg.iloc[i,2])
        string = (start, ':', end)
        id_region.append(' '.join(string))
    
    id_region = pd.DataFrame(list(zip(id_region)), columns = ['region'])
    cov = pd.concat([cvg, id_region], axis = 1)
    
    cov = cov.sort_values(by=['start'])
    cov = cov.drop(['start', 'end'], axis = 1)
    newpath = '{outdir}/coverage'.format(outdir = outdir)
    cov.to_csv(os.path.join(newpath,'{prefix}_coverageresults.tsv'.format(prefix = prefix)), sep = '\t', index = False)

    #Barplot to compare the perbaseavg
    plt.figure(figsize = (30,20))
    sns.set(style = 'white', font_scale = 1.3)
    ax = sns.barplot(x = 'region', y = 'avg_perbase_depth' , hue = 'flag', data = cov, palette = 'RdBu', saturation = 1)
    ax.set_xlabel('Region')
    ax.set_ylabel('Average depth')
    plt.xticks(rotation = 90)
    sns.despine()
    save = ax.get_figure()
    save.savefig(os.path.join(newpath, 'coverage_stats.jpg'))


def cvg_main(mappedlocations, conflictlocations, bamfile, reference, outdir, prefix): 
    """Main function of this script
    
    Parameters
    ----------
    mappedlocations: dataframe
        Coordinates of the mapped regions
    conflictlocations: dataframe
        Coordinates of the conflict regions
    bamfile: str
        The file location of the BAM file
    reference: str
        The file location of the reference FASTA file
    outdir: str
        Output directory
    prefix: str
        Name of genome
        
    """
    newpath = 'coverage'
    os.makedirs(os.path.join(outdir,newpath))
    bar = progressbar.ProgressBar(widgets = ['Running SAMtools: ', progressbar.Bar(), '(', progressbar.ETA(),')'])
    for i in bar(range(1)):
        make_bed(mappedlocations, conflictlocations, reference, outdir, prefix)
        sam(bamfile, outdir, prefix)
        output(outdir, prefix)

