import os
import logging
import argparse
import subprocess
from k2sc import mg2sc

############# Import and set up #############

logging.basicConfig(level=logging.DEBUG)
logging.info("Started the run.")

# Define the CLI defaults:
defaults = {'verbosity': 'debug', 'threads': 2}

# Create the command line interface:
parser = argparse.ArgumentParser(description='Assign metagenomic assignment to single cells')

# required arguments
parser.add_argument('-i', '--input', dest = 'bamfile', help = "Input bam file", required = True)
parser.add_argument('-o', '--outdir', dest = 'outdir', help = "Directory for output", required = True)
parser.add_argument('-db', '--DBpath', dest ='dbfile', help ="Path to kraken database", required = True)

# optional arguments
parser.add_argument('-v', '--verbose', dest = 'verbosity', default = defaults['verbosity'], choices = ['error', 'warning', 'info', 'debug'], help = 'Set logging level (default {verbosity})'.format(**defaults))
parser.add_argument('-n', '--threads', dest = 'threads', default = defaults['threads'], help = "n cores")
parser.add_argument('-prefix', '--prefix', dest = 'prefix', help = "Prefix file output filenames")

# Parse the CLI:
args = parser.parse_args()
args.threads = str(args.threads)

# Set up logging based on the verbosity level set by the command line arguments:
logging.basicConfig(format='%(levelname)s: %(message)s', level=args.verbosity.upper())
    
# extract name prefix from input bam
if args.prefix is not None:
    prefix = args.prefix
else:  
    prefix = os.path.splitext(args.bamfile)[0]

# Generate variables based on input
bamfile_out = os.path.join(args.outdir,prefix + "_unmapped.bam")
fqfile = os.path.join(args.outdir,prefix + "_unmapped.fq")
krakenoutfile = os.path.join(args.outdir,prefix + "_output.kraken")
reportf = os.path.join(args.outdir,prefix + "_krakenreport.txt")

# Make output directories and check that all files exist
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

# Give user feedback about the run
logging.info("Input bamfile: {}".format(args.bamfile))
logging.info("Output krakenfile: {}".format(krakenoutfile))
logging.info("Threads used: {}".format(args.threads))
    
############# Prepare files #############

# Extract unmapped reads from bam
# samtools view -b -f 4 starsoloinput.bam > output_unmapped.bam
cmd1 = "samtools view -b -f 4 " + args.bamfile + " > " + bamfile_out
proc1 = subprocess.Popen(cmd1, shell=True)
proc1.wait()
logging.info("Unmapped reads were extracted and saved to {}".format(bamfile_out))

# Convert to fastq
# bedtools bamtofastq -i output_unmapped.bam -fq output_unmapped.fq
cmd2 = "bedtools bamtofastq -i " + bamfile_out + " -fq " + fqfile
proc2 = subprocess.Popen(cmd2, shell=True)
proc2.wait()
logging.info("FASTQ generated and saved to {}".format(fqfile))

############# Run metagenomic tool on fastq and report summary #############
# run kraken
# kraken2 --threads 2 --db dbpath --report reportfilepath inputfastqpath > krakenoutputfilepath
cmd3 = "kraken2 --threads " + args.threads + " --db " + args.dbfile + " --report " + reportf + " " + fqfile + " > " + krakenoutfile
proc3 = subprocess.Popen(cmd3, shell=True)
proc3.wait()
logging.info("Kraken2 finished running, krakenoutput saved to {}".format(krakenoutfile))

# Split output into single cell level and report sparse matrix (mg2sc.py)
mg2sc(bamfile_out, krakenoutfile, args.dbfile, args.outdir)
logging.info("Sparse matrix with single cell information created")
logging.info("Run finished.")