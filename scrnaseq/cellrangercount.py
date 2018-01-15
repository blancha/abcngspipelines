#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates cellranger count scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=cellranger_count", default="cellranger_count")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files. DEFAULT=../data/FASTQ_files/untrimmed", default="../data/FASTQ_files/untrimmed")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
ppn = config.get("cellrangercount", "ppn")
mem = config.get("cellrangercount", "mem")
vmem = config.get("cellrangercount", "vmem")

genome = config.get("project", "genome")
cellrangerTranscriptome = config.get(genome, "cellrangerTranscriptome")

localcores = config.get("cellrangercount", "localcores")
localmem = config.get("cellrangercount", "localmem")

# Read samples file.
samplesDataFrame = util.readSamplesFile()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Cycle through all the samples and write the cellrangercount scripts.
for index, row in samplesDataFrame.iterrows():
    sample = row["sample"]
    # Create script file.
    scriptName = "cellrangercount_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "cellrangercount")
    script.write("cellranger count " + "\\\n")
    script.write("--localcores=" + localcores + " \\\n")
    script.write("--localmem=" + localmem + " \\\n")
    script.write("--id=" + sample + " \\\n")
    script.write("--fastqs " + os.path.relpath(os.path.join(inputDirectory)) + " \\\n")
    script.write("--sample=" + sample + " \\\n")
    script.write("--transcriptome=" + cellrangerTranscriptome + " \\\n")

    script.write("&> " + scriptName + ".log")    

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
