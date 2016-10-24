#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import os.path
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates dexseqcounts scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="dexseqcounts")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files.", default="../results/tophat")
parser.add_argument("-o", "--outputDirectory", help="Output directory with dexseqcounts results.", default="../results/dexseqcounts")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
toolsFolder = config.get("server", "toolsFolder")

trim = config.getboolean("project", "trim")
stranded = config.getboolean("project", "stranded")
genome = config.get("project", "genome")
genomeFile = config.get(genome, "genomeFile")
dexseq_gtfFile = config.get(genome, "dexseq_gtfFile")
processors = config.get("dexseqcounts", "processors")

# Read samples file.
samplesFile = util.readSamplesFile()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

############################
# dexseqcounts.sh scripts #
############################
for index, row in samplesFile.iterrows():
    sample = row["sample"]
    if "Lane" in samplesFile.columns:
        sample= sample + "_lane_" + str(row["lane"])    
    scriptName = "dexseqcounts_" + sample + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "dexseqcounts")
    script.write("source " + os.path.join(toolsFolder, "python_environments/python2.7/bin/activate"))
    script.write("\n\n")
    script.write("dexseq_count.py" + " \\\n")
    script.write("--paired=yes" + " \\\n")
    if stranded:
        script.write("--stranded=reverse" + " \\\n")
    else:
        script.write("--stranded=no" + " \\\n")
    script.write("--format=bam" + " \\\n")
    script.write("--order=pos" + " \\\n")
    script.write(dexseq_gtfFile + " \\\n")
    script.write(os.path.relpath(os.path.join(inputDirectory, sample, "accepted_hits.bam")) + " \\\n")
    script.write(os.path.relpath(os.path.join(outputDirectory, sample + ".txt")) + " \\\n")
    script.write("&> " + scriptName + ".log")

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
