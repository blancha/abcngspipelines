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
parser = argparse.ArgumentParser(description="Generates Cuffmerge scripts, in the old format.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="cuffmerge")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files.", default="../results/cufflinks")
parser.add_argument("-o", "--outputDirectory", help="Output directory with Tophat results.", default="../results/cuffmerge")
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

genome = config.get("project", "genome")
gtfFile = config.get(genome, "gtfFile")
genomeFile = config.get(genome, "genomeFile")
processors = config.get("cuffmerge", "processors")

# Get samples and conditions
samples = util.getMergedsamples()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Write assembly_GTF_list.txt
assembly_GTF_list = open("assembly_GTF_list.txt", "w")
for sample in samples:
    assembly_GTF_list.write(os.path.abspath(os.path.join("../../results/cufflinks", sample, "transcripts.gtf")) + "\n") 

# Write the cuffmerge script.
scriptName = "cuffmerge.sh"
script = open(scriptName, 'w')
if header:
    util.writeHeader(script, config, "cuffmerge")
script.write("cuffmerge" + " \\\n")
script.write("--ref-gtf " + gtfFile + " \\\n")
script.write("--num-threads " + processors + " \\\n")
script.write("--ref-sequence " + genomeFile + " \\\n") 
script.write("-o " + outputDirectory + " \\\n")
script.write("assembly_GTF_list.txt" + " \\\n")
script.write("&> " + scriptName + ".log")
script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
