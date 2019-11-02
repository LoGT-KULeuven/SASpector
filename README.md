# Integrated Bioinformatics Project 19-20

"Analysis of missing regions in short-read assemblies of bacterial genomes" Emma Verkinderen, Deniz Sinar, Alejandro Correa Rojo.
Advisor: Cédric Loord
Master of Bioinformatics - KU Leuven

## Tool: What the *** is wrong with my Illumina Assembly? (To change)

A bioinformatics tool to extract and analyze unmapped regions from a hybrid genome assembly which were not covered by De novo assembly with Illumina short reads.

## Introduction

What the *** is wrong with my Illumina Assembly? is a tool that compares a Hybrid assembly with a De novo assembly of bacterial genomes by extracting (if apply) unmapped regions from De novo assembly in the hybrid assembly and analyze them to see functional and compositional pattern to explain why these regions are not covered in Illumina sequencing.

The tool takes as global inputs the hybrid assembly as reference genome and short reads assembly as contigs/draft genome, both in FASTA format. What the *** is wrong with my Illumina Assembly? is composed with four modules:

- mapper: mapping of the short reads assembly with the hybrid assembly using progressiveMauve.
- extraction: extraction of the unmapped regions and the regions with ordering problems with the reference genome (Explain this).
- checker: alignment of the unmapped regions and conflicted regions with reference and contigs to detect false positives in the mapping.
- analysis: in progress 

## Getting Started

### Requirements

- Python 3.4 or later
- Java
- Mauve or progressiveMauve (See Bioconda)
- BLAST+
- BioPython
- Pandas

### Installation

in progress

## Usage

### mapper

`mapper3.py` takes the two genomes as input files and one prefix as identification for the output files which will serve as inputs for the other module:

```
python3 mapper3.py reference contigs prefix

usage: What the *** is wrong with my Illumina Assembly? -- Mapper
       [-h] reference contigs prefix

Genome mapper with progressiveMauve for Hybrid and De novo assemblies

positional arguments:
  reference   Hybrid assembly FASTA file as reference genome
  contigs     De novo assembly FASTA file as contigs/draft genome
  prefix      Genome ID

optional arguments:
  -h, --help  show this help message and exit

```

It wraps Mauve and performs an alignment of both genomes with progressiveMauve algorithm. It generates several alignment output files but most importantly the coordinates of the mapped and unmapped regions in the reference genome (the backbone file).

## extract

`extract3`