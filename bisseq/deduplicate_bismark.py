#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import subprocess
import sys
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates deduplicatebismark scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="deduplicatebismark")
parser.add_argument("-e", "--extension", help="Extension to identify SAM files from which duplicates must be removed. DEFAULT=bismark_pe.sam", default="bismark_pe.sam")
parser.add_argument("-i", "--inputDirectory", help="Input directory with SAM files. DEFAULT=../results/bismark", default="../results/bismark")
parser.add_argument("-o", "--outputDirectory", help="Output directory with sorted SAM files. DEFAULT=..results/bismark", default="../results/bismark")
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

# Get samples and conditions
samples = util.getMergedsamples()

# Create scripts directory, if it does not exist yet, and cd to it.
util.makeDirectory(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
util.makeDirectory(outputDirectory)

for sample in samples:
    # Create script
    scriptName = "deduplicatebismark_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "deduplicate_bismark")
    script.write("deduplicate_bismark " + "\\\n")
    script.write("--paired " + "\\\n")
    inputFile = glob.glob(inputDirectory + "/" + sample + "/*" + args.extension)[0]
    script.write(inputFile + " \\\n")
    script.write("&> " + scriptName + ".log")

script.close()

