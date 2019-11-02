from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import pandas as pd
import argparse

def backboneUnmapped(backbonefile):
    backbone = pd.read_table(backbonefile, sep = '\s+')
    unmappedlocations = backbone[(backbone.seq1_leftend == 0) & (backbone.seq1_rightend == 0)]
    unmappedlocations = unmappedlocations[['seq0_leftend', 'seq0_rightend']]
    return unmappedlocations

def backboneConflict(backbonefile):
    backbone = pd.read_table(backbonefile, sep = '\s+')
    conflictlocations = backbone[(backbone.seq0_leftend == 0) & (backbone.seq0_rightend == 0)]
    conflictlocations = conflictlocations[['seq1_leftend', 'seq1_rightend']]
    return conflictlocations

def unmapped(reference, unmappedlocations, prefix):
    unmapped = {}
    idun = []
    for reads in SeqIO.parse(reference, format = 'fasta'):
        for i in range(0, unmappedlocations.shape[0]):
            start = unmappedlocations.iloc[i,0]
            end = unmappedlocations.iloc[i,1]
            unmapped[i] = str(reads.seq[start:end])
    for i in range(0, unmappedlocations.shape[0]):
        start = unmappedlocations.iloc[i,0]
        end = unmappedlocations.iloc[i,1]
        header = (str(prefix), str(start), ':', str(end))
        idun.append(''.join(header))
    for i in range(0, len(unmapped)):
        unmapped[idun[i]] = unmapped.pop(i)
    return unmapped

def conflict(reference, conflictlocations, prefix):
    conflict = {}
    idcon = []
    for reads in SeqIO.parse(reference, format = 'fasta'):
        for i in range(0, conflictlocations.shape[0]):
            start = conflictlocations.iloc[i,0]
            end = conflictlocations.iloc[i,1]
            conflict[i] = str(reads.seq[start:end])
    for i in range(0, conflictlocations.shape[0]):
        start = conflictlocations.iloc[i,0]
        end = conflictlocations.iloc[i,1]
        header = (str(prefix), str(start), ':', str(end))
        idcon.append(''.join(header))
    for i in range(0, len(conflict)):
        conflict[idcon[i]] = conflict.pop(i)
    return conflict

def writeunmapped(unmapped, prefix):
    with open('{prefix}_unmappedregions.fasta'.format(prefix = prefix), 'w+') as fasta:
        for key, value in unmapped.items():
            fasta.write('>' + key + '\n' + value + '\n')

def writeconflict(conflict, prefix):
    with open('{prefix}_conflictcontigs.fasta'.format(prefix = prefix), 'w+') as fasta:
        for key, value in conflict.items():
            fasta.write('>' + key + '\n' + value + '\n')

def main():
    parser = argparse.ArgumentParser(prog = 'What the *** is wrong with my Illumina Assembly? - Extract')
    parser.add_argument('reference', help = 'Hybrid assembly FASTA file as format')
    parser.add_argument('backbone', help = 'Backbone file from progressiveMauve with alignment coordinates')
    parser.add_argument('prefix', help = 'Genome ID')
    args = parser.parse_args()

    unmappedloc = backboneUnmapped(args.backbone)
    conflictloc = backboneConflict(args.backbone)
    unmappeddict = unmapped(args.reference, unmappedloc, args.prefix)
    conflictdict = conflict(args.reference, conflictloc, args.prefix)
    print(unmappeddict)
    writeunmapped(unmappeddict, args.prefix)
    writeconflict(conflictdict, args.prefix)

if __name__ == '__main__':
    main()
