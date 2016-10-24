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
parser = argparse.ArgumentParser(description="Generates meanbedgraphs scripts")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="meanbedgraphs")
parser.add_argument("-i", "--inputDirectory", help="Input directory with bedgraph files. DEFAULT=../results/bedgraph/with_normalisation", default="../results/bedgraph/with_normalisation")
parser.add_argument("-o", "--outputDirectory", help="Output directory where bedgraph files will be written. DEFAULT=../results/bedgraph/with_normalisation", default="../results/bedgraph/with_normalisation")
parser.add_argument("-t", "--stranded", help="Stranded.", choices=["yes", "no", "y", "n"], default="yes")
parser.add_argument("-r", "--sort", help="Should bedgraph files be sorted, before merging and computing mean? DEFAULT=yes.", choices=["yes", "no", "y", "n"], default="yes")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
stranded = args.stranded
if (args.sort == "yes") | (args.sort == "yes"):
    sort = True
else:
    sort = False

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

# Read samples file. Create file first if it doesn't exist.
if not os.path.exists("samples.txt"):
    subprocess.call("samples.py", shell=True)
samplesFile = pandas.read_csv("samples.txt", delim_whitespace=True)

# Get samples
conditions = samplesFile["condition"]
unique_conditions = list(OrderedDict.fromkeys(conditions)) # Remove duplicates.

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

if stranded:
    strands = ["", "_positive", "_negative"]
else:
    strands = [""]

# Cycle through all the conditions and write the meanbedgraphs scripts.
for condition in unique_conditions:
    samples = samplesFile[samplesFile.condition == condition]["sample"].tolist()
    # Create script file.
    for strand in strands:
        scriptName = 'meanbedgraphs_' + condition + strand + '.sh'
        script = open(scriptName, 'w')
        util.writeHeader(script, config, "meanbedgraphs")
        if sort:
            for sample in samples:
                script.write("sort -k 1,1 -k2,2n " + "\\\n")
                script.write(os.path.relpath(os.path.join(inputDirectory, sample + strand + ".bedgraph")) + " \\\n")
                script.write("--output " + os.path.relpath(os.path.join(inputDirectory, sample + strand + ".bedgraph")) + " \\\n")
                script.write("&> " + scriptName + ".log")
                script.write("\n\n")
        script.write("meanbedgraphs_computation.py " + "\\\n")
        script.write("--bedgraphs " + "\\\n")
        for sample in samples[:-1]:
            script.write(os.path.relpath(os.path.join(inputDirectory, sample + strand + ".bedgraph")) + ",\\\n")
        script.write(os.path.relpath(os.path.join(inputDirectory, samples[-1] + strand + ".bedgraph")) + " \\\n")
        script.write("--output " + os.path.relpath(os.path.join(outputDirectory, condition + strand + ".bedgraph")) + " \\\n")
        script.write("&> " + scriptName + ".log")
        script.close()

if (args.submitJobsToQueue == "yes") | (args.submitJobsToQueue == "y"):
    subprocess.call("submitJobs.py", shell=True)
