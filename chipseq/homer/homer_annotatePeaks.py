#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 12/04/2014

import argparse
import glob
import os
import os.path
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates HOMER annotatePeaks scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=homer_annotatePeaks", default="homer_annotatePeaks")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files. DEFAULT=../results/macs_callpeak", default="../results/macs_callpeak")
parser.add_argument("-o", "--outputDirectory", help="Output directory with Tophat results. DEFAULT=../results/homer_annotatePeaks", default="../results/homer_annotatePeaks")
parser.add_argument("-p", "--peaks", help="Narrow or broad peaks. DEFAULT=narrow", choices=["narrow", "broad"], default="narrow")
parser.add_argument("-d", "--directories", help="Samples in subdirectories for each sample. DEFAULT=yes", choices=["yes", "no", "y", "n"], default="no")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
directories = args.directories
peaks = args.peaks

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

genome = config.get('project', 'genome')

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# List comparisons. Normally, should be names of directories in input directory
comparisons = os.listdir(inputDirectory)

# Annotate all the peaks called by macs callpeaks
for comparison in comparisons:
    if (directories == "yes") | (directories == "y"):
        inputDirectory_comparison = os.path.join(inputDirectory, comparison)
    else:   
        inputDirectory_comparison = os.path.join(inputDirectory)
    # Create output directories
    if (directories == "yes") | (directories == "y"):
        outputDirectory_comparison = os.path.join(outputDirectory, comparison)
    else:   
        outputDirectory_comparison = os.path.join(outputDirectory, os.path.splitext(comparison)[0])
    genomeOntologyDirectory = os.path.join(outputDirectory_comparison, "genome_ontology")
    geneOntologyDirectory = os.path.join(outputDirectory_comparison, "gene_ontology")
    if not os.path.exists(outputDirectory_comparison):
        os.makedirs(outputDirectory_comparison)
    if not os.path.exists(genomeOntologyDirectory):  
        os.makedirs(genomeOntologyDirectory)
    if not os.path.exists(geneOntologyDirectory):
        os.makedirs(geneOntologyDirectory)
    # Create script file
    scriptName = 'annotatePeaks_' + comparison + '.sh'
    script = open(scriptName, 'w')
    util.writeHeader(script, config, "annotatePeaks")
    script.write("annotatePeaks.pl " + "\\\n") 
    if (directories == "yes") | (directories == "y"):
        if peaks == "broad":
            script.write(os.path.relpath(os.path.join(inputDirectory_comparison, comparison + "_peaks.broadPeak")) + " \\\n")
        else:
            script.write(os.path.relpath(os.path.join(inputDirectory_comparison, comparison + "_peaks.narrowPeak")) + " \\\n")
    else: 
        script.write(os.path.relpath(os.path.join(inputDirectory_comparison, comparison)) + " \\\n")
    script.write(genome + " -gsize " + genome + " -cons -CpG " + "\\\n")
    script.write("-go " + os.path.relpath(geneOntologyDirectory) + " \\\n")
    script.write("-genomeOntology " + os.path.relpath(genomeOntologyDirectory) + " \\\n")
    script.write("1> " + os.path.relpath(os.path.join(outputDirectory_comparison, "peaks.csv")) + " \\\n")
    script.write("2> " + scriptName + ".log")    
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)

