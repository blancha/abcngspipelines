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
parser = argparse.ArgumentParser(description="Generates Picard tools MarkDuplicates scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=markduplicates", default="markduplicates")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/bismark", default="../results/bismark")
parser.add_argument("-o", "--outputDirectory", help="Output directory with htseqcount results. DEFAULT=../results/bismark", default="../results/bismark")
parser.add_argument("-e", "--extension", help="Extension to identify BAM files in which duplicates must be marked. DEFAULT=.bam", default=".bam")
parser.add_argument("-r", "--removeduplicates", help="Remove duplicates, instead of just marking. DEFAULT=FALSE", default="FALSE")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
extension = args.extension

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

picard_folder = config.get("picard", "folder")
remove_duplicates = config.getboolean("markduplicates", "remove_duplicates")

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
    # Reorder 
    script.write("java -jar " + picard_folder + "/MarkDuplicates.jar " + "\\\n")
    inputFile = glob.glob(inputDirectory + "/" + sample + "/*" + extension)[0]
    script.write("INPUT=" + inputFile + "\\\n")
    #script.write("OUTPUT=" + outputDirectory + "/" + sample + "/" + sample + "_deduplicated.bam " + "\\\n")    
    script.write("OUTPUT=" + outputDirectory + "/" + sample + "_deduplicated.bam " + "\\\n")
    #script.write("METRICS_FILE=" + outputDirectory + "/" + sample + "/" + sample + "_deduplication_metrics.txt " + "\\\n")
    script.write("METRICS_FILE=" + outputDirectory + "/" + sample + "_deduplication_metrics.txt " + "\\\n")
    if remove_duplicates:
        script.write("REMOVE_DUPLICATES=true " + "\\\n")
    else:
        script.write("REMOVE_DUPLICATES=false " + "\\\n")
    script.write("&> " + scriptName + ".log") 

    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
