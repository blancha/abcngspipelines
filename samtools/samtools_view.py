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
parser = argparse.ArgumentParser(description="Generates samtools view scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Directory in which scripts will be placed. Directory will be created, if it doesn't exist yet. DEFAULT=samtools_view", default="samtools_view")
parser.add_argument("-e", "--extension", help="Extension to identify BAM files to view. DEFAULT=.bam", default=".bam")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=..results/tophat", default="../results/tophat")
parser.add_argument("-o", "--outputDirectory", help="Output directory with viewed BAM files. DEFAULT=..results/tophat", default="../results/tophat")
parser.add_argument("-r", "--region", help="Region to view")
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

extension = args.extension
region = args.region

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

processors = config.get("samtools", "processors")

# Get samples and conditions
samples = util.getMergedsamples()
samples = sorted(samples, reverse=True) # To always put wt first, put the list in reverse alpabetical order

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

for sample in samples:
    # Create script
    scriptName = "samtools_view_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "samtoolsIndex")
    script.write("samtools view -h " + "\\\n")
    inputFile = glob.glob(inputDirectory + "/" + sample + "/*" + extension)[0]
    script.write(inputFile + " \\\n")
    script.write(region + " \\\n")
    #script.write("1> " + outputDirectory + "/" + sample + "/" + sample + "_filtered_" + extension + " \\\n")
    script.write("1> " + outputDirectory + "/" + sample + "/" + sample + "_ERCC_only.bam" + " \\\n")
    script.write("2> " + scriptName + ".log")

script.close()

