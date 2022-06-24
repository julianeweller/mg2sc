# mg2sc
Analyzing metagenomics of single cell RNA-seq

The aim of the project was developing a method that extracts unmapped reads and uses
metagenomic tools to classify their taxonomy on a single-cell level. This information is
quantified for each transcript and cell, resulting in a count matrix with cell ID by transcript count for each organism.

# Usage
Developed on python 3.8.

**Required python packages**
- pysam v0.16.0.1
- scipy v1.6.2
- regex v2021.4.4   

**Required command line packages**
- kraken2 (https://github.com/DerrickWood/kraken2, conda install kraken2)
- samtools (https://github.com/samtools)
- bedtools v2.30.0 (https://bedtools.readthedocs.io)

**Setting up kraken2**

In additional to these packages, you'll need to setup a reference database for kraken2. Please see the kraken2 manual for how to do this (https://github.com/DerrickWood/kraken2/wiki/Manual). Alternatively, you can download pre-built Kraken 2 database (e.g. Standard from 12/2/2020, 36GB) from https://benlangmead.github.io/aws-indexes/k2. 

Please ensure that you have enough memory available for reading in the kraken2 database, e.g. for the example database, more than 40 GB of memory is recommended. If you have less memory available, you could go with a smaller reference database.

**Execution**
Both files, scMeG-kraken.py and k2sc.py are required for the full workflow and need to be in the same folder.

command:\
python scMeG-kraken.py --input [bamfile, e.g. starsolo/Aligned.sortedByCoord.out.bam] --outdir [output directory] / --DBpath [path to kraken database] --threads [#, e.g. 8] --prefix [prefered file prefix] --verbosity [error/warning/info/debug]

successfull run:
- "Sparse matrix with single cell information created"
- "Run finished."

**Interpretation of results**

The pipeline produces a cellranger2-style formatted output folder, which can be easily imported for downstream analysis. The observation space is the complete list of barcodes with reads identified by the pipeline, and the feature space is complete list of organisms (and other elements of the kraken2 hierarchy) that were identified.

It can be useful to collapse this count matrix to more general categories, such as family or genus. You can use `src/collapse_taxonomy.py` for that, with a quick demo shown in [a notebook](demo/collapse_taxonomy_demo.ipynb).

# Background

The first step of the tool workflow is the file preparation. From the user input, that
most importantly includes the read alignment file and the
corresponding reference database, the unaligned reads are extracted, and the data is
converted into a FASTQ file.

The FASTQ file serves as input for the metagenomic analysis. The output, sequencing read IDs with assigned taxonomy, can be combined with the cell barcode and transcript UMI from the alignment file with unmapped reads to assign each cell and each transcript a taxonomy. If a transcript has different taxonomies assigned, the taxonomy ID with the highest read count is chosen. 

Finally, the results can be exported as a sparse matrix, which facilitates the integration with the differential gene expression data and cell type annotation into AnnData objects. For Kraken2, the read is assigned to the lowest common ancestor that it maps to. Therefore, the sparse matrix contains the lowest mapped taxonomic level.

![scMeG-kraken](https://user-images.githubusercontent.com/46549848/121743588-7a275e80-caf9-11eb-8d6c-82fb1c217cff.png)

