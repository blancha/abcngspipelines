#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description='Generate mergesrafastqfiles scripts.')
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=mergesrafastqfiles", default="mergesrafastqfiles")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files to merge. DEFAULT=.", default=".")
parser.add_argument("-o", "--outputDirectory", help="Output directory with merged FASTQ files. DEFAULT=.", default=".")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
scriptsDirectory = os.path.abspath(args.scriptsDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

# Get samples
samplesFile = util.readSamplesFile(inputDirectory)

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Keep only the unique sample names
uniquesamplesFile = samplesFile.drop_duplicates(subset="sample")

# Write the scripts
for index, row in uniquesamplesFile.iterrows():
    uniquesample = row["sample"]
    samples = samplesFile[samplesFile.sample == uniquesample]
    # Create output directory, if it does not exist yet..
    if not os.path.exists(os.path.join(inputDirectory, uniquesample)):
        os.makedirs(os.path.join(inputDirectory, uniquesample))
    # Create script
    scriptName = "mergesrafastqfiles_" + uniquesample + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "mergesrafastqfiles")
    script.write("cat" + " \\\n")
    script.write(os.path.join(outputDirectory, uniquesample, uniquesample + ".bam") + " \\\n")
    for index2, row2 in samples.iterrows():
        script.write(os.path.join(inputDirectory, row2["sample"] + "_lane_" + str(row2["Lane"]), row2["sample"] + "_lane_" + str(row2["Lane"]) + "_recal_reads.bam") + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
