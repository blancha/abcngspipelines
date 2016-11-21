#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

# Remove NAs from MATS_output
# for k in *; do echo $k; cd ${k}/MATS_output; for i in *; do j=$(basename $i .txt)_no_na.txt; echo $j; grep -v NA $i > $j; done; cd ../..; done

import argparse
from collections import OrderedDict
import glob
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates rmats2sashimiplots2sashimiplots scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="rmats2sashimiplots")
parser.add_argument("-b", "--bamDirectory", help="Input directory with BAM files.", default="../results/star")
parser.add_argument("-m", "--matsDirectory", help="Input directory with MATS events files.", default="../results/mats")
parser.add_argument("-o", "--outputDirectory", help="Output directory with rmats2sashimiplots2sashimiplots results.", default="../results/rmats2sashimiplots")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
bamDirectory = os.path.abspath(args.bamDirectory)
matsDirectory = os.path.abspath(args.matsDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Check if the bamDirectory exists, and is a directory.
util.checkInputDirectory(bamDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
toolsFolder = config.get("server", "toolsFolder")

if config.has_option("rmats2sashimiplots", "comparisons"):
    comparisons = config.get("rmats2sashimiplots", "comparisons")
else:
    comparisons = config.get("mats", "comparisons")

# Read samples file. Create file first if it doesn't exist.
if not os.path.exists("samples.txt"):
    subprocess.call("samples", shell=True)
samplesFile = pandas.read_csv("samples.txt", delim_whitespace=True)

# Get samples
samples = list(samplesFile["sample"])
conditions = list(samplesFile["group"])
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

events = ["A3SS", "A5SS", "MXE", "RI", "SE"]

countTypes = ["JunctionCountOnly", "ReadsOnTargetAndJunctionCounts"]

# Write the rmats2sashimiplots script.
for comparison in comparisons:
    comparison = comparison.replace(" ", "_").lower()
    # Create comparisons subdirectory in scripts directory, if it does not exist yet, and cd to it.
    if not os.path.exists(comparison):
        os.mkdir(comparison)
    os.chdir(comparison)  
    # Create comparison subdirectory in output directory, if it does not exist yet.
    if not os.path.exists(os.path.join(outputDirectory, comparison)):
        os.mkdir(os.path.join(outputDirectory, comparison))
    condition1 = comparison.split("_vs_")[0].strip()
    samples1 =  list(samplesFile[samplesFile["group"] == condition1]["sample"])
    condition2 = comparison.split("_vs_")[1].strip()
    samples2 =  list(samplesFile[samplesFile["group"] == condition2]["sample"])
    # Create script for each event
    for event in events:
         # Create event subdirectory in comparison directory, if it does not exist yet.
        if not os.path.exists(os.path.join(outputDirectory, comparison, event)):
            os.mkdir(os.path.join(outputDirectory, comparison, event))       
            # Create script for each count type.
        for countType in countTypes:
             # Create countType subdirectory in event directory, if it does not exist yet.
            if not os.path.exists(os.path.join(outputDirectory, comparison, event, countType.lower())):
                os.mkdir(os.path.join(outputDirectory, comparison, event, countType.lower()))
            scriptName = "rmats2sashimiplots_" + comparison + "_" + event + "_" + countType.lower() + ".sh"
            script = open(scriptName, 'w')
            if header:
                util.writeHeader(script, config, "rmats2sashimiplots")
            script.write("# Deactivate Python 3 virtual environment, and activate Python 2 virtual environment" + "\n")    
            script.write("source " + os.path.join(toolsFolder, "python_environments/python2.7/bin/activate") + "\n")
            script.write("\n")
            script.write("rmats2sashimiplot" + " \\\n")
            script.write("-b1 ")
            for sample in samples1[:-1]:
                script.write(os.path.relpath(os.path.join(bamDirectory, sample, sample + ".bam")) + ",")
            script.write(os.path.relpath(os.path.join(bamDirectory, samples1[-1], samples1[-1] + ".bam")) + " \\\n")
            script.write("-b2 ")
            for sample in samples2[:-1]:
                script.write(os.path.relpath(os.path.join(bamDirectory, sample, sample + ".bam")) + ",")
            script.write(os.path.relpath(os.path.join(bamDirectory, samples2[-1], samples2[-1] + ".bam")) + " \\\n")
            script.write("-t SE" + " \\\n")
            script.write("-e " + os.path.relpath(os.path.join(matsDirectory, comparison, "MATS_output", event + ".MATS." + countType + "_no_na.txt")) + " \\\n")
            script.write("-l1 " + condition1 + " \\\n")
            script.write("-l2 " + condition2 +" \\\n")
            script.write("-o " + os.path.relpath(os.path.join(outputDirectory, comparison, event, countType.lower())) + " \\\n")
            script.write("&> " + scriptName + ".log")
            script.close()
    os.chdir("..")

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs", shell=True)
