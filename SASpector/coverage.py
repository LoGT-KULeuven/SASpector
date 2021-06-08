# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 21:58:29 2019

@author: emmav, 0mician
"""
import logging

import pandas as pd
import csv
from Bio import SeqIO
import subprocess
import os
import seaborn as sns
import matplotlib.pyplot as plt

""" coverage

This script calculates the average coverage of the mapped and unmapped regions. Using SAMtools bedcov, it calculates
the coverage for both regions and generates a summary table including the coordinates for each region, total read base count
and coverage defined as average per base depth. Additionally, it generates a boxplot of the average per base depth for each mapped and unmapped region.
The input file is a provided alignment (BAM file) of the short reads against the reference genome.

"""

# bwa mem {inputRef} {input.R1} {input.R2} | samtools view -b -o -- | bwa index-a bwtsw {inputREF}

def make_bed(mappedlocations, conflictlocations, reference, outdir, prefix):
    """ Generates bed files for the mapped locations, unmapped locations and conflict locations.
        For coverage analysis the unmapped regions are not filtered and flanks are not included.

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
    logging.info("Bed file generation")
    for seq in SeqIO.parse(reference, "fasta"):
        ID = seq.id.split(' ')[0]

    #Include only the filtered (> 100 bp) unmapped regions from summary file
    unmap = []
    unmapsum = '{outdir}/{prefix}_unmapsummary.tsv'.format(outdir = outdir, prefix = prefix)
    regions = pd.read_csv(unmapsum, sep = '\t')
    regions = regions['Region'].values.tolist()
    for region in regions:
        start = region.split('_')[1].split(':')[0]
        end = region.split('_')[1].split(':')[1]
        dictunmap = {'id': ID, 'start': start, 'end': end}
        unmap.append(dictunmap)

    #Mapped regions: from dataframes generated in summary script
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
    """Wraps SAMtools to sort and index the BAM file, then determine the coverage for both mapped and unmapped regions with samtools bedcov.
       It calculates the total read base count (sum of per base depth) per region. The result tsv files are saved in new subdirectory 'coverage'.

    Parameters
    ----------
    bamfile: str
        The file location of the BAM file
    outdir: str
        Output directory
    prefix: str
        Name of genome

    """
    logging.info("Samtools analysis")
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
    """Writes the coverage statistics for each mapped and unmapped region to a result tsv file, and generates a boxplot (jpg) of the coverage for each region.
       Coverage is defined as average per base depth over the region.

    Parameters
    ----------
    outdir: str
        Output directory
    prefix: str
        Name of genome
    """
    logging.info("Saving coverage analysis files")
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

    cvg = cvg.sort_values(by=['start'])
    newpath = '{outdir}/coverage'.format(outdir = outdir)
    cvg.to_csv(os.path.join(newpath,'{prefix}_coverageresults.tsv'.format(prefix = prefix)), sep = '\t', index = False)

    #Boxplot to compare the perbaseavg
    plt.figure()
    sns.set(style = 'white', font_scale = 1.3)
    ax = sns.boxplot(x = 'flag', y = 'avg_perbase_depth', data = cvg, fliersize=2.5)
    ax.set_xlabel('')
    ax.set_ylabel('Average depth')
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles, ["Mapped", "Unmapped"])
    sns.despine()
    save = ax.get_figure()
    save.savefig(os.path.join(newpath, 'coverage_boxplots.jpg'))
    plt.clf()


def cvg_main(mappedlocations, conflictlocations, bamfile, reference, outdir, prefix):
    """Main function of this script

    Parameters
    ----------
    mappedlocations: dataframe
        Coordinates of the mapped regions in the reference genome
    conflictlocations: dataframe
        Coordinates of the conflict regions in the reference genome
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
    logging.info("Running coverage analysis (bed file generation, samfile analysis, saving files)")
    make_bed(mappedlocations, conflictlocations, reference, outdir, prefix)
    sam(bamfile, outdir, prefix)
    output(outdir, prefix)
    logging.info("Coverage analysis completed!")
