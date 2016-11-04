#!/usr/bin/env python3

# Version 1.2
# Author Alexis Blanchet-Cohen
# Last modified: 27/10/2016

import argparse
import glob
import os
import os.path
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description='Generates FastQC scripts.')
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=fastqc", default="fastqc")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files. DEFAULT=../data/FASTQ_files/untrimmed", default="../data/FASTQ_files/untrimmed")
parser.add_argument("-o", "--outputDirectory", help="Output directory with FastQC results. DEFAULT=../results/fastqc", default="../results/fastqc")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Read samples file.
samplesDataFrame = util.readSamplesFile()
samples = samplesDataFrame["sample"].tolist()

# Create script and output directories, if they do not exist yet.
util.makeDirectory(outputDirectory, recursive=True)
util.makeDirectory(scriptsDirectory, recursive=True)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

# Change to scripts directory
os.chdir(scriptsDirectory)

files = os.listdir(inputDirectory)

# Cycle through all the files and write the fastqc scripts.
for file in files:
    if not os.path.isfile(file):
	continue
    scriptName = "fastqc_" + os.path.basename(file) + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "fastqc")
    script.write("fastqc " + "\\\n")
    script.write("--outdir " + os.path.relpath(outputDirectory)  + " \\\n")
    script.write(os.path.relpath(os.path.join(file)) + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)

