#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
from collections import OrderedDict
import glob
import os
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates Cufflinks scripts, in the old format.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="cufflinks")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files.", default="../results/tophat")
parser.add_argument("-o", "--outputDirectory", help="Output directory with Tophat results.", default="../results/cufflinks")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

stranded = config.getboolean("project", "stranded")
genome = config.get("project", "genome")
genomeFile = config.get(genome, "genomeFile")
gtfFile = config.get(genome, "gtfFile")
processors = config.get("cufflinks", "processors")

# Get samples and conditions
samples = util.getMergedsamples()
#samples = sorted(samples, reverse=True) # To always put wt first, put the list in reverse alpabetical order
conditions = util.getMergedconditions()
#conditions = sorted(conditions, reverse=True) # To always put wt first, put the list in reverse alpabetical order
unique_conditions = list(OrderedDict.fromkeys(conditions)) # Remove duplicates.

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Write the cufflinks script.
for sample in samples:
    if not os.path.exists(outputDirectory + "/" + sample):
        os.mkdir(outputDirectory + "/" + sample)
    scriptName = "cufflinks_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "cufflinks")
    script.write("cufflinks" + " \\\n")
    script.write("--num-threads " + processors + " \\\n")
    script.write("-o " + outputDirectory + "/" + sample + " \\\n")
    script.write(inputDirectory + "/" + samples[0] + "/accepted_hits.bam" + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.close()
    if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
        subprocess.call("submitJobs.py", shell=True)
