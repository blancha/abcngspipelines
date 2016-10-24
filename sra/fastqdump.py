#!/usr/bin/env python3


# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates scripts to fastq-dump SRA files.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=fastqdump", default="fastqdump")
parser.add_argument("-a", "--samplesFile", help="Input file with names of SRA runs. DEFAULT=SraRunTable.txt", default="SraRunTable.txt")
parser.add_argument("-i", "--inputDirectory", help="Input directory with SRA files. DEFAULT=../data/SRA_files", default="../data/SRA_files")
parser.add_argument("-o", "--outputDirectory", help="Output directory with FASTQ files. DEFAULT=.", default="../data/FASTQ_files/untrimmed")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
samplesFile = os.path.abspath(args.samplesFile)
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Check if the samplesFile exists, and is a file.
if not(os.path.exists(samplesFile) and os.path.isfile(samplesFile)):
    exit(samplesFile +  " does not exist or is not a file.")

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

# Read input file.
samplesFile = pandas.read_csv(samplesFile, sep="\t")

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Cycle through all the samples and write the star scripts.
for index, row in samplesFile.iterrows():
    run = row["Run_s"]
    # Create script file.
    scriptName = "fastqdump_" + run + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "getsra")
    script.write("fastq-dump" + " \\\n")
    script.write("--gzip" + " \\\n")
    script.write("--split-files" + " \\\n")
    script.write("--outdir " + os.path.relpath(outputDirectory) + " \\\n")
    script.write(os.path.relpath(os.path.join(inputDirectory, run + ".sra")) + " \\\n")

    script.write("&> " + scriptName + ".log")    

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
