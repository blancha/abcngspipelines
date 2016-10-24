#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
from collections import OrderedDict
import glob
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates MATS scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="mats")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files.", default="../results/star")
parser.add_argument("-o", "--outputDirectory", help="Output directory with MATS results.", default="../results/mats")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
toolsFolder = config.get("server", "toolsFolder")

stranded = config.getboolean("project", "stranded")
genome = config.get("project", "genome")
gtfFile = config.get(genome, "gtfFile")
comparisons = config.get("mats", "comparisons")

# Read samples file. Create file first if it doesn't exist.
if not os.path.exists("samples.txt"):
    subprocess.call("samples.py", shell=True)
samplesFile = pandas.read_csv("samples.txt", delim_whitespace=True)

# Get samples
samples = list(samplesFile["sample"])
conditions = list(samplesFile["condition"])
unique_conditions = list(OrderedDict.fromkeys(conditions)) # Remove duplicates

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

def processComparisons(comparisons):
    comparisons = comparisons.split(",")
    comparisons = [str.strip(x) for x in comparisons]
    return(comparisons)        

comparisons = processComparisons(comparisons)

# Write the mats script.
for comparison in comparisons:
    if not os.path.exists(os.path.relpath(os.path.join(outputDirectory, comparison.replace(" ", "_").lower()))):
        os.mkdir(os.path.relpath(os.path.join(outputDirectory, comparison.replace(" ", "_").lower())))
    condition1 = comparison.split("vs")[0].strip()
    samples1 =  list(samplesFile[samplesFile["condition"] == condition1]["sample"])
    condition2 = comparison.split("vs")[1].strip()
    samples2 =  list(samplesFile[samplesFile["condition"] == condition2]["sample"])
    scriptName = "mats_" + comparison.replace(" ", "_").lower() + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "mats")

    script.write("# Deactivate Python 3 virtual environment, and activate Python 2 virtual environment" + "\n")    
    script.write("source " + os.path.join(toolsFolder, "python_environments/python2.7/bin/activate") + "\n")
    script.write("\n")
    script.write("RNASeq-MATS.py" + " \\\n")
    script.write("-b1 ")
    for sample in samples1[:-1]:
        script.write(os.path.relpath(os.path.join(inputDirectory, sample, sample + ".bam")) + ",")
    script.write(os.path.relpath(os.path.join(inputDirectory, samples1[-1], samples1[-1] + ".bam")) + " \\\n")
    script.write("-b2 ")
    for sample in samples2[:-1]:
        script.write(os.path.relpath(os.path.join(inputDirectory, sample, sample + ".bam")) + ",")
    script.write(os.path.relpath(os.path.join(inputDirectory, samples2[-1], samples2[-1] + ".bam")) + " \\\n")
    script.write("-gtf " + gtfFile + " \\\n")
    script.write("-o " + os.path.relpath(os.path.join(outputDirectory, comparison.replace(" ", "_").lower())) + " \\\n")
    script.write("-t paired" + " \\\n")
    script.write("-len 50" + " \\\n")
    script.write("-r1 150,150,150" + " \\\n")
    script.write("-r2 150,150,150" + " \\\n")
    script.write("-sd1 75,75,75" + " \\\n")
    script.write("-sd2 75,75,75" + " \\\n")
    script.write("-analysis U" + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
