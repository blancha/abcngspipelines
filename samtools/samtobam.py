#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import os.path
import subprocess
import sys
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates samtobam (samtools view) scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="samtobam")
parser.add_argument("-e", "--extension", help="Extension to identify SAM files to convert to BAM files. DEFAULT=.sam", default=".bam")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/bowtie", default="../results/bowtie")
parser.add_argument("-o", "--outputDirectory", help="Output directory with sorted BAM files. DEFAULT=..results/bowtie", default="../results/bowtie")
#parser.add_argument("-u", "--subDirectories", help="samples in subdirectories. DEFAULT=..yes", choices=["yes", "no", "y", "n"], default="yes")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the scripts directory, cd to the scripts directory
util.cdMainScriptsDirectory()

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

# Read samples file.
samplesFile = util.readSamplesFile()

# Create scripts directory, if it does not exist yet, and cd to it.
util.makeDirectory(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
util.makeDirectory(outputDirectory)

for index, row in samplesFile.iterrows():
    sample = row["sample"]
    # Create script
    scriptName = "samtobam_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "samtoolsIndex")
    script.write("samtools view -hb " + os.path.relpath(os.path.join(outputDirectory, sample, sample + ".sam")) + " \\\n")
    script.write("1> " + os.path.relpath(os.path.join(outputDirectory + "/" + sample + "/" + sample + ".bam")) + " \\\n")
    script.write("2> " + scriptName + ".log")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)


