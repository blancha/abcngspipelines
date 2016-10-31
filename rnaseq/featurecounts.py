#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates featurecounts scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="featurecounts")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files.", default="../results/star")
parser.add_argument("-o", "--outputDirectory", help="Output directory with featurecounts results.", default="../results/featurecounts")
parser.add_argument("-e", "--extension", help="Extension of BAM files on which to run featureCounts", default=".bam")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
extension = args.extension

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

trim = config.getboolean("project", "trim")
stranded = config.getboolean("project", "stranded")
genome = config.get("project", "genome")
genomeFile = config.get(genome, "genomeFile")
gtfFile = config.get(genome, "gtfFile")
T = config.get("featurecounts", "T")
p = config.getboolean("featurecounts", "p")
B = config.getboolean("featurecounts", "B")
C = config.getboolean("featurecounts", "C")
s = config.getboolean("featurecounts", "s")
M = config.getboolean("featurecounts", "M")
O = config.getboolean("featurecounts", "O")

# Read samples file.
samplesFile = util.readSamplesFile()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

############################
# featurecounts.sh scripts #
############################
for index, row in samplesFile.iterrows():
    sample = row["sample"]
    if "lane" in samplesFile.columns:
        sample += "_" + row["lane"]   
    scriptName = "featurecounts_" + sample + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "featurecounts")
    script.write("featureCounts \\\n")
    script.write("-T " + T + " \\\n")
    if p:
        script.write("-p" + " \\\n")
    if B:
        script.write("-B" + " \\\n")
    if C:
        script.write("-C" + " \\\n")
    if s:
        script.write("-s 2" + " \\\n")
    if M:
        script.write("-M" + " \\\n")
    if O:
        script.write("-O" + " \\\n")

    script.write("-a " + gtfFile + " \\\n")
    script.write("-o " + os.path.relpath(os.path.join(outputDirectory, sample + ".txt")) + " \\\n")
    script.write(os.path.relpath(os.path.join(inputDirectory, sample, sample + extension)) + " \\\n")
    script.write("&> " + scriptName + ".log") 

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
