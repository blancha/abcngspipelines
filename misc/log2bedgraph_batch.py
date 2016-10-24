#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 21/03/2014

import argparse
import configparser
import glob
import os
import pandas
import subprocess
import sys
import util

parser = argparse.ArgumentParser(description='Generate scripts to normalise bedgraph files by scaling factor.')
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=log2bedgraph", default="log2bedgraph")
parser.add_argument("-i", "--inputDirectory", help="Input directory with bedgraph files. DEFAULT=../results/bedgraph/normalised", default="../results/bedgraph/normalised")
parser.add_argument("-o", "--outputDirectory", help="Output directory with normalised bedgraph files. DEFAULT=../results/bedgraph/log2", default="../results/bedgraph/log2")
parser.add_argument("-a", "--sampleSubDirectories", help="BAM files located in sample subdirectories in inputDirectory. DEFAULT=yes", choices=["yes", "no", "y", "n"], default="yes")
parser.add_argument("-e", "--extension", help="Extension identifying bedgraph files. DEFAULT=.bedgraph", default=".bedgraph") 
parser.add_argument("-t", "--stranded", help="Stranded bedgraph files. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
if (args.stranded == "yes") | (args.stranded == "y"):
    stranded = True
else:
    stranded = False

# Read samples file. Create file first if it doesn't exist.
if not os.path.exists("samples.txt"):
    subprocess.call("samples.py", shell=True)
samplesFile = pandas.read_csv("samples.txt", sep="\t")

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(inputDirectory)

# Create script and output directories, if they do not exist yet.
util.makeDirectory(outputDirectory)
util.makeDirectory(scriptsDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
genome = config.get("project", "genome")

# cd to scripts directory
os.chdir(scriptsDirectory)

# Cycle through all the samples and write the tophat scripts.
for index, row in samplesFile.iterrows():
    sample = row["sample"]
    # Create script file.
    scriptName = 'log2bedgraph_' + sample + '.sh'
    script = open(scriptName, 'w')
    util.writeHeader(script, config, "log2bedgraph_batch")
    script.write("log2bedgraph.py " + "\\\n")
    script.write("--before_log2 " + inputDirectory + "/" + sample + ".bedgraph " + "\\\n")
    script.write("--after_log2 " + outputDirectory + "/" + sample + ".bedgraph " + "\\\n")
    script.write("&> " + scriptName + ".log")
    script.close()
    if stranded: 
        # Create positive strand script file.
        scriptName = 'log2bedgraph_' + sample + '_positive.sh'
        script = open(scriptName, 'w')
        util.writeHeader(script, config, "log2bedgraph_batch")
        script.write("log2bedgraph.py " + "\\\n")
        script.write("--before_log2 " + inputDirectory + "/" + sample + "_positive.bedgraph " + "\\\n")
        script.write("--after_log2 " + outputDirectory + "/" + sample + "_positive.bedgraph " + "\\\n")
        script.write("&> " + scriptName + ".log")
        script.close()
        # Create negative strand script file.
        scriptName = 'log2bedgraph_' + sample + '_negative.sh'
        script = open(scriptName, 'w')
        util.writeHeader(script, config, "log2bedgraph_batch")
        script.write("log2bedgraph.py " + "\\\n")
        script.write("--before_log2 " + inputDirectory + "/" + sample + "_negative.bedgraph " + "\\\n")
        script.write("--after_log2 " + outputDirectory + "/" + sample + "_negative.bedgraph " + "\\\n")
        script.write("&> " + scriptName + ".log")
        script.close()
if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
