#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 13/04/2014

import argparse
import configparser
import glob
import os
import os.path
import subprocess
import sys
import util

parser = argparse.ArgumentParser(description="Separates BAM files by strand.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=separatebamsbystrand", default="separatebamsbystrand")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/star", default="../results/star")
parser.add_argument("-o", "--outputDirectory", help="Output directory with bedgraph files. DEFAULT=../results/star", default="../results/star")
parser.add_argument("-e", "--extension", help="Extension to identify BAM files to sort. DEFAULT=.bam", default=".bam")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Read samples files
samples = util.getSamples()

# Process command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
extension = args.extension

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(inputDirectory)

# Create script and output directories, if they do not exist yet.
util.makeDirectory(outputDirectory)
util.makeDirectory(scriptsDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

# cd to scripts directory
os.chdir(scriptsDirectory)

# Function to check BAM file names.
# Create a symbolic link to accepted_hits.bam named sampleName.bam. Useful for TopHat output, which is named accepted_hits.bam, by default
#def createSymbolicLinks():
#    for sample in samples:
#        command = "ln -fs " + os.path.join(inputDirectory, sample, "accepted_hits.bam") + " " + os.path.join(inputDirectory, sample, sample + ".bam")
#        print(command)
#        subprocess.call(command, shell=True)

#if (args.symbolicLinks.lower() == "yes") | (args.symbolicLinks.lower() == "y"):
#    createSymbolicLinks()

strands = ["positive", "negative"]

# Write script for each sample and each strand
for sample in samples:
    for strand in strands:
        # Create script file.
        scriptName = "separatebamsbystrand_" + sample + "_" + strand + ".sh"
        script = open(scriptName, "w")
        if header:
            util.writeHeader(script, config, "separatebambystrand")
        # BAM to bedgraph
        script.write("separatebambystrand.py" + " \\\n")
        script.write("--input_bam " + os.path.relpath(os.path.join(inputDirectory, sample, sample + ".bam")) + " \\\n")
        script.write("--strand " + strand + " \\\n")
        script.write("&> " + scriptName + ".log")

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
