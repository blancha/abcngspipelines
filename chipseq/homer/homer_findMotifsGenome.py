#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 12/04/2014

import argparse
import glob
import os
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates HOMER findMotifs scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=homer_findMotifs", default="homer_findMotifs")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files. DEFAULT=../results/macs_callpeak", default="../results/macs_callpeak")
parser.add_argument("-o", "--outputDirectory", help="Output directory with Tophat results. DEFAULT=../results/homer_findMotifs", default="../results/homer_findMotifs")
parser.add_argument("-p", "--peaks", help="Narrow or broad peaks. DEFAULT=narrow", choices=["narrow", "broad"], default="narrow")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
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

# List comparison. Normally, should be names of directories in input directory
comparisons = os.listdir(inputDirectory)

# Annotate all the peaks called by macs callpeaks
for comparison in comparisons:
    inputDirectory_comparison = os.path.join(inputDirectory, comparison)
    # Create output directories
    outputDirectory_comparison = os.path.join(outputDirectory, comparison)
    if not os.path.exists(outputDirectory_comparison):
        os.makedirs(outputDirectory_comparison)
    # Create script file
    scriptName = 'findMotifs_' + comparison + '.sh'
    script = open(scriptName, 'w')
    util.writeHeader(script, config, "findMotifs")
    script.write("findMotifsGenome.pl " + "\\\n")  
    if peaks == "narrow":
        script.write(os.path.join(inputDirectory_comparison, comparison + "_peaks.narrowPeak") + " \\\n")
    if peaks == "broad":
        script.write(os.path.join(inputDirectory_comparison, comparison + "_peaks.broadPeak") + " \\\n")
    script.write(genome + " \\\n")
    script.write(os.path.join(outputDirectory_comparison) + " \\\n")
    script.write("-size 150 " + "\\\n")
    script.write("&> " + scriptName + ".log")    
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
