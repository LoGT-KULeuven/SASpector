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

# bwa mem {inputRef} {input.R1} {input.R2} | samtools view -b -o -- | bwa index-a bwtsw {inputREF}


#make bedfiles for mapped and unmapped regions
def make_bed(mappedlocations, conflictlocations, reference, outdir, prefix):
    for seq in SeqIO.parse(reference, "fasta"):
        ID = seq.id.split(' ')[0]

    #only the filtered unmapped regions from summary file
    unmap = []
    unmapsum = '{outdir}/{prefix}_unmapsummary.tsv'.format(outdir = outdir, prefix = prefix)
    regions = pd.read_table(unmapsum)['Region'].values.tolist()
    for region in regions:
        start = region.split('_')[1].split(':')[0]
        end = region.split('_')[1].split(':')[1]
        dictunmap = {'id': ID, 'start': start, 'end': end}
        unmap.append(dictunmap)

    #mapped regions: from df generated in summary script
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
    sort = 'samtools sort -o {outdir}/coverage/{prefix}.sorted.bam {bam} &&'.format(outdir = outdir, prefix = prefix, bam = bamfile)
    index = 'samtools index {outdir}/coverage/{prefix}.sorted.bam &&'.format(outdir = outdir, prefix = prefix)
    bedcovu = 'samtools bedcov {outdir}/coverage/{prefix}_unmappedregions.bed {outdir}/coverage/{prefix}.sorted.bam > {outdir}/coverage/{prefix}_unmapcvg.tsv &&'.format(outdir = outdir, prefix = prefix)
    bedcovm = 'samtools bedcov {outdir}/coverage/{prefix}_mappedregions.bed {outdir}/coverage/{prefix}.sorted.bam > {outdir}/coverage/{prefix}_mapcvg.tsv'.format(outdir = outdir, prefix = prefix)
   
    pipe1 = subprocess.Popen(sort, stdout = subprocess.PIPE)
    pipe2 = subprocess.Popen(index, stdin = pipe1.stdout, stdout = subprocess.PIPE)
    pipe3 = subprocess.Popen(bedcovu, stdin = pipe2.stdout, stdout = subprocess.PIPE)
    pipe4 = subprocess.Popen(bedcovm, stdin = pipe3.stdout, stdout = subprocess.PIPE)
    output = pipe4.communicate()
    pipe1.stdout.close()
    print(output.decode())



def output(prefix, outdir):
    bedcovu = '{outdir}/coverage/{prefix}_unmapcvg.tsv'.format(prefix = prefix, outdir = outdir)
    bedcovm = '{outdir}/coverage/{prefix}_mapcvg.tsv'.format(prefix = prefix, outdir = outdir)
    cvgu = pd.read_table(bedcovu, header=None)
    cvgu['perbaseavg'] = cvgu[3]/(cvgu[2]-cvgu[1])
    #compare with mapped regions
    cvgm = pd.read_table(bedcovm, header=None)
    cvgm['perbaseavg'] = cvgm[3]/(cvgm[2]-cvgm[1])

    cvgu['flag']='unmapped'
    cvgm['flag']='mapped'
    cvg=cvgu.append(cvgm)
    cvg.columns=['seq','start', 'end', 'total_perbase_depth', 'avg_perbase_depth', 'flag']
    cvg.sort_values(by=['start'])
    
    #barplot to compare the perbaseavg


def cvg_main(mappedlocations, conflictlocations, bamfile, reference, outdir, prefix):
    newpath = 'coverage'
    os.makedirs(os.path.join(outdir,newpath))

    bar = progressbar.ProgressBar(widgets = ['Running SAMtools: ', progressbar.Bar(), '(', progressbar.ETA(),')'])
    for i in bar(range(1)):
        sam(bamfile, outdir, prefix)
        make_bed(mappedlocations, conflictlocations, reference, outdir, prefix)
        output(prefix, outdir)
