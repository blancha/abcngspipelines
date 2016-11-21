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
parser = argparse.ArgumentParser(description="Generates Cuffdiff scripts, in the old format.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="cuffdiff_with_cuffquant")
parser.add_argument("-i", "--inputDirectory", help="Input directory with CXB files.", default="../results/cuffquant")
parser.add_argument("-o", "--outputDirectory", help="Output directory with Cuffdiff results.", default="../results/cuffdiff")
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

stranded = config.getboolean("project", "stranded")
genome = config.get("project", "genome")
genomeFile = config.get(genome, "genomeFile")
gtfFile = config.get(genome, "gtfFile")
processors = config.get("cuffdiff_with_cuffquant", "processors")

# Read samples file. Create file first if it doesn't exist.
if not os.path.exists("samples.txt"):
    subprocess.call("samples.py", shell=True)
samplesFile = pandas.read_csv("samples.txt", delim_whitespace=True)

# Get samples
samples = samplesFile["sample"]
conditions = samplesFile["group"]
unique_conditions = list(OrderedDict.fromkeys(conditions)) # Remove duplicates.

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Write the cuffdiff script.
scriptName = "cuffdiff_with_cuffquant.sh"
script = open(scriptName, 'w')
if header:
    util.writeHeader(script, config, "cuffdiff")
script.write("cuffdiff" + " \\\n")
script.write("--labels ")
script.write(",".join(unique_conditions) + " \\\n")
script.write("-p " + processors + " \\\n")
script.write("--no-effective-length-correction " + "\\\n")
if stranded:
    script.write("--library-type fr-firststrand" + " \\\n")
script.write("-u -b " + genomeFile + " \\\n") 
script.write("-o " + os.path.relpath(os.path.join(outputDirectory)) + " \\\n")
script.write(gtfFile + " \\\n")
script.write(os.path.relpath(os.path.join(inputDirectory, samples[0], "abundances.cxb")))
previous_condition = conditions[0]
for sample,condition in zip(samples[1:],conditions[1:]):
    if condition == previous_condition:
        script.write(",\\\n")
    else:
        script.write(" \\\n")
    script.write(os.path.relpath(os.path.join(inputDirectory, sample, "abundances.cxb")))
    previous_condition = condition
script.write(" \\\n")
script.write("&> " + scriptName + ".log")
script.close()

# Copy the CummeRbund script.
command = "cp -r " + os.path.join(config.get("server", "toolsFolder"), "abcngspipelines/rnaseq/cummerbund") +  " " + os.path.dirname(scriptsDirectory)
print(command)
subprocess.call(command, shell=True)

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
