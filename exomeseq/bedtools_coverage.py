#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates bedtools coverage scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=bedtools_coverage", default="bedtools_coverage")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files. DEFAULT=../results/bwa", default="../results/bwa")
parser.add_argument("-o", "--outputDirectory", help="Output directory with bedtools coverage results. DEFAULT=../results/bedtools_coverage", default="../results/bedtools_coverage")
parser.add_argument("-a", help="One or more BAM/BED/GFF/VCF file(s) “B”. Use “stdin” if passing B with a UNIX pipe. NEW!!!: -b may be followed with multiple databases and/or wildcard (*) character(s). DEFAULT=/gs/project/feb-684-aa/BIF/VIN/VIN-BIF-P1/data/SeqCapEZ_Exome_v3.0_Design_Annotation_files/120430_HG19_ExomeV3_UTR_EZ_HX1_ensembl_sorted_filtered.bed", 
        default="/gs/project/feb-684-aa/BIF/VIN/VIN-BIF-P1/data/SeqCapEZ_Exome_v3.0_Design_Annotation_files/120430_HG19_ExomeV3_UTR_EZ_HX1_ensembl_sorted_filtered.bed")
parser.add_argument("-c", "--chromSizes", help="Tab delimited file with chromosome sizes and lengths. DEFAULT=/gs/project/feb-684-aa/BIF/genomes/Homo_sapiens/Broad/human_g1k_v37_chrom.sizes", default="/gs/project/feb-684-aa/BIF/genomes/Homo_sapiens/Broad/human_g1k_v37_chrom.sizes")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
a = os.path.abspath(args.a)
chromSizes = os.path.abspath(args.chromSizes)

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

# Read samples file.
samplesFile = util.readsamplesFile()
samples = samplesFile["sample"].tolist()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Cycle through all the samples and write the bedtools_coverage scripts.
for sample in samples:
    # Create script file.
    scriptName = "bedtools_coverage_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "bedtools_coverage")
    script.write("bedtools coverage" + " \\\n")
    script.write("-a " + os.path.relpath(a) + " \\\n")
    script.write("-b " + os.path.relpath(os.path.join(inputDirectory, sample, sample + ".bam")) + " \\\n")
    script.write("-g " + os.path.relpath(chromSizes) + " \\\n")
    script.write("-hist" + " \\\n")
    script.write("-sorted" + " \\\n")
    script.write("1> " + os.path.relpath(os.path.join(outputDirectory, sample + ".txt")) + " \\\n")
    script.write("2> " + scriptName + ".log")    

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
