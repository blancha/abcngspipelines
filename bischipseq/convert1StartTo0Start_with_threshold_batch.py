#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 15/06/2014

import argparse
import glob
import os
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description='Generate scripts to convert bedgraph files from one-based start to zero-based start.')
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="convert1StartTo0Start_with_threshold")
parser.add_argument("-i", "--inputDirectory", help="Input directory with bedgraph files.", default="../bedgraph/methylation_counts_sorted/")
parser.add_argument("-o", "--outputDirectory", help="Output directory with sorted bedgraph files.", default="../bedgraph/methylation_counts_sorted_0_start/")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

samples = util.getMergedsamples()

# Read configuration files.
config = util.readConfigurationFiles()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.mkdir(outputDirectory)

for file in os.listdir(inputDirectory):
    file = os.path.splitext(file)[0]
    # Create script file.
    scriptName = 'convert1StartTo0Start_with_threshold_' + file + '.sh'
    script = open(scriptName, 'w')
    util.writeHeader(script, config, "convert1StartTo0Start_with_threshold")
    script.write("convert1StartTo0Start_with_threshold.py " + "\\\n")
    script.write("--one_start_bedgraph " + inputDirectory +  "/" + file + ".bedgraph " + "\\\n")
    script.write("--zero_start_bedgraph_with_threshold " + outputDirectory + "/" + file + ".bedgraph")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
