# mg2sc
Analyzing metagenomics of single cell RNA-seq

The aim of the project was developing a method that extracts unmapped reads and uses
metagenomic tools to classify their taxonomy on a single-cell level. This information is
quantified for each transcript and cell, resulting in a count matrix with cell ID by transcript count for each organism.

# Usage

## Required python packages
- pysam v0.16.0.1
- scipy v1.6.2
- regex v2021.4.4   

## Required command line packages
- kraken2 (https://github.com/DerrickWood/kraken2, conda install kraken2)
- samtools (https://github.com/samtools)
- bedtools v2.30.0 (https://bedtools.readthedocs.io)

## Setting up kraken2

In additional to these packages, you'll need to setup a reference database for kraken2. Please see the kraken2 manual for how to do this (https://github.com/DerrickWood/kraken2/wiki/Manual). Alternatively, you can download pre-built Kraken 2 database (e.g. Standard from 12/2/2020, 36GB) from https://benlangmead.github.io/aws-indexes/k2

## Execution
command:\
python scMeG-kraken.py --input [bamfile, e.g. starsolo/Aligned.sortedByCoord.out.bam] --outdir [output directory] / --DBpath [path to kraken database] --threads [#, e.g. 8] --prefix [prefered file prefix] --verbosity [error/warning/info/debug]

successfull run:
- "Sparse matrix with single cell information created"
- "Run finished."

## Interpretation of results

# Background
The first step of the tool workflow is the file preparation. From the user input, that
most importantly includes the read alignment file and the
corresponding reference database, the unaligned reads are extracted, and the data is
converted into a FASTQ file.

The FASTQ file serves as input for the metagenomic analysis. The output, sequencing read IDs with assigned taxonomy, can be combined with the cell barcode and transcript UMI from the alignment file with unmapped reads to assign each cell and each transcript a taxonomy. If a transcript has different taxonomies assigned, the taxonomy ID with the highest read count is chosen. 

Finally, the results can be exported as a sparse matrix, which facilitates the integration with the differential gene expression data and cell type annotation into AnnData objects. For Kraken2, the read is assigned to the lowest common ancestor that it maps to. Therefore, the sparse matrix contains the lowest mapped taxonomic level.

![image](https://user-images.githubusercontent.com/46549848/121742917-692a1d80-caf8-11eb-98e2-8a709a3512a9.png)
