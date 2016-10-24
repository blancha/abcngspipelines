#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates bwa mem scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="bwa")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files.", default="../data/FASTQ_files/untrimmed/")
parser.add_argument("-o", "--outputDirectory", help="Output directory with bwa results.", default="../results/bwa/")
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
genome = config.get("project", "genome")
genomeFile = config.get(genome, "genomeFile")
bwaIndex = config.get(genome, "bwaIndex")
threads = config.get("bwamem", "threads")

# Read samples file
samplesFile = util.readsamplesFile()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directories, if they do not exist yet..
if not os.path.exists(outputDirectory): 
    os.makedirs(outputDirectory)

# Cycle through all the samples and write the bwa scripts.
for index, row in samplesFile.iterrows():
    sample = row["sample"]
    if "lane" in samplesFile.columns:
        sample= sample + "_lane_" + str(row["lane"])
# Create output directories
    if not os.path.exists(outputDirectory + "/" + sample):
        os.mkdir(outputDirectory +"/" + sample)
    file_R1 = row["file_r1"]
    file_R2 = row["file_r2"]
    # Create script file.
    scriptName = 'bwa_' + sample + '.sh'
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "bwamem")
    script.write("bwa mem -M" + " \\\n")
    script.write("-t " + threads + " \\\n")
    script.write("-R '@RG\\tID:" + sample + "\\tSM:" + row["sample"] + "\\tPL:Illumina\\tLB:lib1\\tPU:unit1'" + " \\\n");
    script.write(bwaIndex + " \\\n")
    script.write(inputDirectory + "/" + file_R1 + " \\\n")
    script.write(inputDirectory + "/" + file_R2 + " \\\n")
    script.write("1> " + outputDirectory + "/" +  sample + "/" + sample + ".sam " + "\\\n")
    script.write("2> " + scriptName + ".log") 

    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
