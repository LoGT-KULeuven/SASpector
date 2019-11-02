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

### extract

`extract3.py` takes the hybrid assembly genome (reference) as input FASTA file, the backbone file generated in `mapper3.py` and one prefix as identification for the output files.

```
python3 extract3.py reference backbone prefix

usage: What the *** is wrong with my Illumina Assembly? - Extract
       [-h] reference backbone prefix

Extraction of the missing regions with their locations in the reference genome
for futher analysis. If some regions were not correctly mapped, they will be
extracted as "conflictcontigs."

positional arguments:
  reference   Hybrid assembly FASTA file as format
  backbone    Backbone file from progressiveMauve with alignment coordinates
  prefix      Genome ID

optional arguments:
  -h, --help  show this help message and exit

```
It will parse the backbone file and extract the sequences which are not covered in the short reads assembly from the reference genome and create at multi-fasta file with the coordinates. Additionally, if some regions were not mapped correctly due to gaps or indels, it will also extract them as 'conflict contigs' into multi-fasta file with their locations for further analysis.

### check

`check3.py` is composed by two commands: **blast2Unmap** which is going to align the unmapped regions with the contigs in order to detect is those regions are repeats already included in the contigs from the short reads assembly. Similar, **blast2Conf** will detect if regions from the short reads assembly that where not correctly mapped are already included in the reference genome.

```
python3 check3.py {blast2Unmap,blast2Conf} prefix

usage: What the *** is wrong with my Illumina Assembly? - Checker
       [-h] {blast2Unmap,blast2Conf} ... prefix

Multiple alignment of unmapped regions and conflict contigs to detect false
positives in the genome mapping

positional arguments:
  {blast2Unmap,blast2Conf}
    blast2Unmap         Alignment of unmapped regions
    blast2Conf          Alignment of conflict contigs
  prefix

optional arguments:
  -h, --help            show this help message and exit

```
Both commands will generate a BLAST output alignment results in tab-delimited format.

### analysis