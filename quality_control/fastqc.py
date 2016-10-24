#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 16/05/2014

import argparse
import glob
import os
import os.path
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description='Generates FastQC scripts.')
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=fastqc", default="fastqc")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FastQC files. DEFAULT=../data/FASTQ_files/untrimmed", default="../data/FASTQ_files/untrimmed")
parser.add_argument("-o", "--outputDirectory", help="Output directory with FastQC results. DEFAULT=../results/fastqc", default="../results/fastqc")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
#util.cdMainScriptsDirectory()

# Process command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(inputDirectory)

# Create script and output directories, if they do not exist yet.
util.makeDirectory(outputDirectory, recursive=True)
util.makeDirectory(scriptsDirectory, recursive=True)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

# Change to scripts directory
os.chdir(scriptsDirectory)

# Store all the FASTQ_files filenames.
files = glob.glob(inputDirectory + "/*fastq") + glob.glob(inputDirectory + "/*fastq.gz")

# Write script
for file in files:                  
    scriptName = "fastqc_" + os.path.basename(file) + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "fastqc")
    script.write("fastqc " + "\\\n")
    script.write("--outdir " + os.path.relpath(outputDirectory)  + " \\\n")
    script.write(file + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)

