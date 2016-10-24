#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 17/06/2014

import argparse
from collections import OrderedDict
import os
import subprocess
import util

parser = argparse.ArgumentParser(description='Converts bedgraph files to bigWig files.')
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=bedgraphtobigwig", default="bedgraphtobigwig")
parser.add_argument("-i", "--inputDirectory", help="Input directory with bedgraph files. default=../results/bedgraph/with_normalisation/", default="../results/bedgraph/with_normalisation/")
parser.add_argument("-o", "--outputDirectory", help="Output directory with bigWig files. default=../results/bigwig/with_normalisation", default="../results/bigwig/with_normalisation")
parser.add_argument("-t", "--stranded", help="Stranded files. DEFAULT=yes", choices=["yes", "no", "y", "n"], default="yes")
parser.add_argument("-r", "--sort", help="Sort bedgraph files. If the files are not sorted, they will be sorted with bedtoolSort first. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
stranded = args.stranded
sort = args.sort

# Read samples file.
samplesFile = util.readSamplesFile()

# Get samples
samples = samplesFile["sample"]
conditions = samplesFile["condition"]
unique_conditions = list(OrderedDict.fromkeys(conditions)) # Remove duplicates.

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

genome = config.get("project", "genome")
genomeFolder = config.get(genome, "genomeFolder")
institute = config.get(genome, "institute")

# Create script directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

if stranded:
    strands = ["", "_positive", "_negative"]
else:
    strands = [""]

for sample in samples:
    for strand in strands:
        # Create script file.
        scriptName = 'bedgraphtobigwig_' + sample + strand + '.sh'
        script = open(scriptName, 'w')
        if header:
            util.writeHeader(script, config, "bedgraphtobigwig")
        if sort == "yes":
            script.write("bedSort " + "\\\n") 
            script.write(os.path.relpath(os.path.join(inputDirectory, sample + strand + ".bedgraph")) + " \\\n")
            script.write(os.path.relpath(os.path.join(inputDirectory, sample + strand + ".bedgraph")) + " \\\n")
            script.write("&> " + scriptName + ".log")
            script.write("\n\n")
        script.write("bedGraphToBigWig " + "\\\n")
        script.write(os.path.relpath(os.path.join(inputDirectory, sample + strand + ".bedgraph")) + " \\\n")
        if institute == "Ensembl":
            script.write(os.path.join(genomeFolder, "Annotation", "Genes", "UCSCChromInfo.txt") + " \\\n")
        else: 
            script.write(os.path.join(genomeFolder, "Annotation", "Genes", "ChromInfo.txt") + " \\\n")
        script.write(os.path.relpath(os.path.join(outputDirectory, sample + strand + ".bw")) + " \\\n")
        if sort == "yes":
            script.write("&>> " + scriptName + ".log")
        else:
            script.write("&> " + scriptName + ".log")
        script.close()

for condition in unique_conditions:
    for strand in strands:
        # Create script file.
        scriptName = 'bedgraphtobigwig_' + condition + strand + '.sh'
        script = open(scriptName, 'w')
        if header:
            util.writeHeader(script, config, "bedgraphtobigwig")
        if sorted == "yes":
            script.write("bedSort " + "\\\n") 
            script.write(os.path.relpath(os.path.join(inputDirectory, condition + strand + ".bedgraph")) + " \\\n")
            script.write(os.path.relpath(os.path.join(inputDirectory, condition + strand + ".bedgraph")) + " \\\n")
            script.write("&> " + scriptName + ".log")
            script.write("\n\n")
        script.write("bedGraphToBigWig " + "\\\n")
        script.write(os.path.relpath(os.path.join(inputDirectory, condition + strand + ".bedgraph")) + " \\\n")
        if institute == "Ensembl":
            script.write(os.path.join(genomeFolder, "Annotation", "Genes", "UCSCChromInfo.txt") + " \\\n")
        else: 
            script.write(os.path.join(genomeFolder, "Annotation", "Genes", "ChromInfo.txt") + " \\\n")
        script.write(os.path.relpath(os.path.join(outputDirectory, condition + strand + ".bw")) + " \\\n")
        if sorted == "no":
            script.write("&>> " + scriptName + ".log")
        else:
            script.write("&> " + scriptName + ".log")

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)

