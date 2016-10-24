#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 13/04/2014

import argparse
import configparser
import glob
import os
import subprocess
import sys
import util

parser = argparse.ArgumentParser(description='Converts BAM files to bigwig files.')
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=bamtobigwig", default="bamtobigwig")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/star", default="../results/star")
parser.add_argument("-e", "--bedgraphDirectoryEnsembl", help="Output directory with intermediary Ensembl bedgraph files. DEFAULT=../results/bedgraph/ensembl", default="../results/bedgraph/ensembl")
parser.add_argument("-c", "--bedgraphDirectoryUCSC", help="Output directory with intermediary UCSC bedgraph files. DEFAULT=../results/bedgraph/ucsc", default="../results/bedgraph/ucsc")
parser.add_argument("-o", "--outputDirectory", help="Output directory with bigwig files. DEFAULT=../results/bigwig", default="../results/bigwig/without_normalisation")
parser.add_argument("-t", "--stranded", help="Stranded files.", choices=["yes", "no", "y", "n"], default="yes")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
inputDirectory = os.path.abspath(args.inputDirectory)
bedgraphDirectoryEnsembl = os.path.abspath(args.bedgraphDirectoryEnsembl)
bedgraphDirectoryUCSC = os.path.abspath(args.bedgraphDirectoryUCSC)
outputDirectory = os.path.abspath(args.outputDirectory)
stranded = args.stranded

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(inputDirectory)

# Get samples
samples = util.getSamples()

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
genome = config.get("project", "genome")
genomeFolder = config.get(genome, "genomeFolder")
institute = config.get(genome, "institute")

# Create script and output directories, if they do not exist yet.
util.makeDirectory(bedgraphDirectoryUCSC, recursive=True)
if institute == "Ensembl":
    util.makeDirectory(bedgraphDirectoryEnsembl, recursive=True)
util.makeDirectory(outputDirectory, recursive=True)
util.makeDirectory(scriptsDirectory)

# cd to scripts directory
os.chdir(scriptsDirectory)

if stranded:
    strands = ["", "_positive", "_negative"]
else:
    strands = [""]

# Write script for each file
for sample in samples:
    for strand in strands:
        # Create script file.
        scriptName = 'bamtobigwig_' + sample + strand + '.sh'
        script = open(scriptName, 'w')
        if header:
            util.writeHeader(script, config, "bamtobigwig")
        # BAM to bedgraph
        script.write("bedtools genomecov" + " \\\n")
        # Only if type == rnaseq ???
        script.write("-split" + " \\\n")
        script.write("-bg" + " \\\n")
        script.write("-g " + os.path.join(genomeFolder, "Sequence", "WholeGenomeFasta", "chrom.sizes") + " \\\n")
        script.write("-ibam " + os.path.relpath(os.path.join(inputDirectory, sample, sample + strand + ".bam")) + " \\\n")
        if institute == "Ensembl":
            script.write("1> " + os.path.relpath(os.path.join(bedgraphDirectoryEnsembl, sample + strand + ".bedgraph")) + " \\\n")
        else:
            script.write("1> " + os.path.relpath(os.path.join(bedgraphDirectoryUCSC, sample + strand + ".bedgraph")) + " \\\n")
        script.write("2> " + scriptName + ".log")
        script.write("\n\n")
        # Ensembl to UCSC conversion, only if institute is Ensembl
        if institute == "Ensembl":
            script.write("ensemblbedgraphtoucscbedgraph.py" + " \\\n")
            script.write("--ensembl_bedgraph " + os.path.relpath(os.path.join(bedgraphDirectoryEnsembl, sample + strand + ".bedgraph")) + " \\\n")
            script.write("--ucsc_bedgraph " + os.path.relpath(os.path.join(bedgraphDirectoryUCSC, sample + strand + ".bedgraph")) + " \\\n")
            script.write("--dictionary " + os.path.join(genomeFolder, "Annotation", "Genes", "ucsc_to_ensembl_chromosomes.dict") + " \\\n")
            script.write("&>> " + scriptName + ".log")
            script.write("\n\n")
        # bedgraph to bigwig
        script.write("bedGraphToBigWig " + "\\\n")
        script.write(os.path.relpath(os.path.join(bedgraphDirectoryUCSC, sample + strand + ".bedgraph")) + " \\\n")
        if institute == "Ensembl":
            script.write(os.path.join(genomeFolder, "Annotation", "Genes", "UCSCChromInfo.txt") + " \\\n")
        else:
            script.write(os.path.join(genomeFolder, "Annotation", "Genes", "ChromInfo.txt") + " \\\n")
        script.write(os.path.relpath(os.path.join(outputDirectory, sample + strand + ".bw")) + " \\\n")
        script.write("&>> " + scriptName + ".log")
        script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
