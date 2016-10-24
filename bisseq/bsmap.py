#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 06/03/2014

import argparse
import os
import subprocess
import sys
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description='Generate BSMAP scripts.')
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="bsmap")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files.", default="../data/FASTQ_files/trimmed/")
parser.add_argument("-o", "--outputDirectory", help="Output directory with bsmap results.", default="../results/bsmap/")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Read samples file
samplesFile = open("samples.txt")
samplesFile = samplesFile.readlines()[1:]

# Read configuration files.
config = util.readConfigurationFiles()

genome = config.get("project", "genome")
genomeFile = config.get(genome, "genomeFile")
processors = config.get("bsmap", "processors")

# Check if samples are in more than 1 lane
multipleLanes = util.multipleLanes()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.mkdir(outputDirectory)

# Cycle through all the samples and write the bsmap scripts.
for line in samplesFile:
    sample = line.split()[2]
    if multipleLanes:
        lane = line.split()[3]
        sample= sample + "_" + lane
    # Create output directory for the sample.
    if not os.path.exists(outputDirectory + "/" + sample):
        os.mkdir(outputDirectory + "/" +  sample)    
    file_R1 = line.split()[0]
    file_R2 = line.split()[1]
    # Create script file.
    scriptName = 'bsmap_' + sample + '.sh'
    script = open(scriptName, 'w')
    util.writeHeader(script, config, "bsmap")
    script.write("bsmap -r 0 -s 16 -n 1" + " \\\n")
    script.write("-a " + inputDirectory + "/" + file_R1 + " \\\n")
    script.write("-b " + inputDirectory + "/" + file_R2 + " \\\n")
    script.write("-d " + genomeFile + " \\\n")
    script.write("-p " + processors + " \\\n") 
    script.write("-o " + outputDirectory + "/" + sample + ".sam" + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
