#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import subprocess
import sys
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates samtools sort scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Directory in which scripts will be placed. Directory will be created, if it doesn't exist yet. DEFAULT=samtools_sort_by_coordinates", default="samtools_sort_by_coordinates")
parser.add_argument("-b", "--by", help="Sort by coordinates or read name. DEFAULT=coordinates", choices=["coordinates", "name"], default="coordinates")
parser.add_argument("-e", "--extension", help="Extension to identify BAM files to sort. DEFAULT=.bam", default=".bam")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/star", default="../results/star")
parser.add_argument("-o", "--outputDirectory", help="Output directory with sorted BAM files. DEFAULT=../results/star", default="../results/star")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the scripts directory, cd to the scripts directory
util.cdMainScriptsDirectory()

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
extension_before_sorting = args.extension

by = args.by
if by == "name":
    scriptsDirectory.replace("coordinates", "name")

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

processors = config.get("samtools", "processors")

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

if by == "name":
    extension_after_sorting = "_sorted_by_name"
else:
    extension_after_sorting = "_sorted_by_coordinates"

for sample in samples:
    # Create script
    scriptName = "samtools_sort_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "samtoolsIndex")
    script.write("samtools sort " + "\\\n")
    if by == "name":
        script.write("-n " + "\\\n")
    script.write("-o " + os.path.relpath(os.path.join(outputDirectory, sample, sample + extension_after_sorting + ".bam")) +  " \\\n")
    bamFile = glob.glob(os.path.join(inputDirectory, sample, "*" + extension_before_sorting))[0]
    script.write(os.path.relpath(bamFile) + " \\\n")
    script.write("&> " + scriptName + ".log")

    script.write("\n\n")

    script.write("samtools index " + os.path.relpath(os.path.join(inputDirectory, sample,sample + extension_after_sorting + ".bam"))  +  " \\\n")
    script.write("&>> " + scriptName + ".log")

script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
