# Integrated Bioinformatics Project 19-20

"Analysis of missing regions in short-read assemblies of bacterial genomes" Emma Verkinderen, Deniz Sinar, Alejandro Correa Rojo.

Advisor: Cédric Lood

Master of Bioinformatics - KU Leuven, Belgium

## Tool: SASpector (Short-read Assembly inSpector)

A bioinformatics tool to extract and analyze missing regions of short-read assemblies by mapping the contigs to a reference genome.

## Introduction

SASpector is a tool that compares a short-read assembly with a reference bacterial genome (for example obtained via hybrid assembly) by extracting missing (unmapped) regions from the reference and analyzing them to see functional and compositional pattern. The aim of the analysis is to explain why these regions are missed by the short-read assembly and if important parts of the genome are missed when a resolved genome is lacking.

The tool takes as global inputs the reference genome and a short-read assembly as contigs/draft genome, both in FASTA format. This repository contains a main script `SASpector.py` to obtain missing regions and several other python scripts that can be run separately by the user for evaluation and analysis.

- `mapper.py`: mapping of the short-read assembly against the reference assembly using progressiveMauve. 
- `summary.py`: extraction of the mapped, unmapped (missing) and conflict regions to fasta files. Also creates summary statistics for the unmapped regions and reference which are written to separate csv files.
- `check.py`: running BLAST alignment of the missing regions to the contigs and the conflict regions to the reference to detect false positives in the mapping. 

Optionally, some scripts for analysis of the unmapped regions can be run: 
- `quastunmap.py`: summary statistics using QUAST.
- `kmer.py`: k-mer analysis using KAT.


## Getting Started

### Requirements

- Python 3.4 or later
- Java
- Seaborn
- Mauve or progressiveMauve (See Bioconda)
- BLAST+
- BioPython
- Pandas 
- KAT (optional)
- QUAST (optional)

### Installation

`git clone` this repository

## Usage

For basic functionalities, run `SASpector.py` with Python3. This wraps the mapping and extraction steps (`mapper.py` and `summary.py`)
```
usage: SASpector [-h] [-p PREFIX] [-dir OUTDIR] [-f [FLANKING]] reference contigs

positional arguments:
  reference             Reference genome FASTA file
  contigs               Short-read assembly (draft genome) FASTA file

optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        Genome ID
  -dir OUTDIR, --outdir OUTDIR
                        Output directory
  -f [FLANKING], --flanking [FLANKING]
                        Add flanking regions to the extracted sequences
                        [Default = 0 bp]
                        
  ```

First, Mauve performs an alignment of both genomes with the progressiveMauve algorithm. It will generate a subdirectory prefix.alignment with several output files but most importantly the backbone file with coordinates of the mapped and unmapped regions in the reference genome. 

Afterwards, this script will parse the backbone file and extract the sequences that are not covered in the short-read assembly from the reference genome. They are written to a multi-fasta file with the prefix and coordinates in the headers, which is done equally for the mapped and conflict regions (regions that didn't map correctly due to gaps or indels). 

Finally, two tab-delimited summary files are generated in a subdirectory called summary. One for the reference, with the amount of gapped and ungapped regions, the fraction of the reference genome that they represent, the GC content and the length. The other one for the unmapped regions, with for each region the GC content and length and then for each amino acid the occurence frequency averaged over all 6 reading frames.

### Checker - in progress (blast settings tuning)

`check.py` consists of two functions: blast2Unmap and blast2Conf, both using blastn. **blast2Unmap** aligns the unmapped regions to the contigs in order to check if those regions are included in the contigs from the short-read assembly. If yes, they could be repeats. Similarly, **blast2Conf** will detect if the conflict regions are included in the reference genome.

```
python3 check.py {blast2Unmap,blast2Conf} prefix

usage: SASpector - Checker
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
Both commands will generate a BLAST result file in tab-delimited format.

### Analysis

The following python scripts are provided for analysis of missing regions. 
- `quastunmap.py`
- `kmer.py`
- ... (also add some explanations)

They are not (yet) implemented as command-line tools, anyone is free to use them at will. For the IBP project, they are used in a Snakemake file for automation and high-throughput goals.

### References
- Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. Journal of Molecular Biology, 215(3), 403–410.
- Darling, A. C. E. (2004). Mauve: Multiple Alignment of Conserved Genomic Sequence With Rearrangements. Genome Research, 14(7), 1394–1403. 
- Mapleson, D., Garcia Accinelli, G., Kettleborough, G., Wright, J., & Clavijo, B. J. (2016). KAT: a K-mer analysis toolkit to quality control NGS datasets and genome assemblies. Bioinformatics, 33(4).
- Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: Quality assessment tool for genome assemblies. Bioinformatics, 29(8), 1072–1075. 
- Sth else?
