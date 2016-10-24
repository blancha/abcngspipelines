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
parser = argparse.ArgumentParser(description="Generates Picard tools CollectInsertMetrics scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=picard_collectinsertmetrics", default="picard_collectinsertmetrics")
parser.add_argument("-e", "--extension", help="Extension to identify BAM files on which to collect metrics. DEFAULT=.bam", default=".bam")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/star", default="../results/star")
parser.add_argument("-o", "--outputDirectory", help="Output directory. DEFAULT=../results/picard_collectinsertmetrics", default="../results/picard_collectinsertmetrics")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
extension=args.extension

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
picardFolder = config.get("picard", "folder")

# Read samples file.
samplesFile = util.readSamplesFile()
samples = samplesFile["sample"].tolist()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Cycle through all the samples and write the picard_collectinsertmetrics scripts.
for sample in samples:
    # Create script file.
    scriptName = "picard_collectinsertmetrics_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "picard_collectinsertmetrics")
    script.write("java -Xmx2700m" + " \\\n")
    script.write("-jar " + os.path.join(picardFolder, "picard.jar") + " \\\n")
    script.write("CollectInsertSizeMetrics "+ " \\\n")
    script.write("I=" + os.path.relpath(os.path.join(inputDirectory, sample, sample + extension)) + " \\\n")
    script.write("O=" + os.path.relpath(os.path.join(outputDirectory, sample + "insert_size_metrics.txt")) + " \\\n")
    script.write("H=" + os.path.relpath(os.path.join(outputDirectory, sample + "insert_size_histogram.pdf")) + " \\\n")
    script.write("&> " + scriptName + ".log")    

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
