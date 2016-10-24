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
parser = argparse.ArgumentParser(description="Generates Picard tools MarkDuplicates scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=markduplicates", default="markduplicates")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/bismark", default="../results/bwa")
parser.add_argument("-o", "--outputDirectory", help="Output directory with htseqcount results. DEFAULT=../results/bwa", default="../results/bwa")
parser.add_argument("-r", "--removeduplicates", help="Remove duplicates, instead of just marking. DEFAULT=FALSE", default="FALSE")
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

picard_folder = config.get("picard", "folder")
xmx = config.get("picard", "xmx")

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
    if "Lane" in samplesFile.columns:
        sample= sample + "_lane_" + str(row["Lane"])
    scriptName =  "markduplicates_" + sample + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "markduplicates")
    # Sort 
    script.write("java -Xmx" + xmx + " \\\n")
    script.write("-jar " + picard_folder + "/SortSam.jar" + " \\\n")
    script.write("INPUT=" + os.path.join(inputDirectory, sample, sample + ".sam") + " \\\n")
    script.write("OUTPUT=" + os.path.join(outputDirectory, sample, sample + "_sorted.bam") + " \\\n")    
    script.write("SORT_ORDER=coordinate" + " \\\n")     
    script.write("&> " + scriptName + ".log")
    script.write("\n\n")    
    # Mark duplicates
    script.write("java -Xmx" + xmx + " \\\n")
    script.write("-jar " + picard_folder + "/MarkDuplicates.jar" + " \\\n")
    script.write("INPUT=" + os.path.join(inputDirectory, sample, sample + "_sorted.bam")  + " \\\n")
    script.write("OUTPUT=" + os.path.join(outputDirectory, sample, sample + "_deduplicated.bam") + " \\\n")
    script.write("METRICS_FILE=" + os.path.join(outputDirectory, sample, sample + "_deduplication_metrics.txt") + " \\\n")
    script.write("&>> " + scriptName + ".log")
    script.write("\n\n") 
    # Index BAM file
    script.write("java -Xmx" + xmx + " \\\n")
    script.write("-jar " + picard_folder + "/BuildBamIndex.jar" + " \\\n")
    script.write("INPUT=" + os.path.join(outputDirectory, sample, sample + "_deduplicated.bam") + " \\\n")
    script.write("&>> " + scriptName + ".log")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
