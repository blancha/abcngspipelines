#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 21/05/2014

import argparse
import glob
import os
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description='Sorts bedgraph files, with the following command: sort -k1,1 -k2,2n.')
parser.add_argument("-c", "--scriptsDirectory", help="Scripts directory.", default="sortbedgraph")
parser.add_argument("-i", "--inputDirectory", help="Input directory with bedgraph files.", default="../bedgraph/methylation_counts/")
parser.add_argument("-o", "--outputDirectory", help="Output directory with sorted bedgraph files.", default="../bedgraph/methylation_counts_sorted/")
parser.add_argument("-s", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

samples = util.getMergedSamples()

# Read configuration files.
config = util.readConfigurationFiles()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.mkdir(outputDirectory)

files = os.listdir(inputDirectory)
for file in files:
    # Create script file.
    file = os.path.splitext(file)[0]
    scriptName = 'sortbedgraph_' + file + '.sh'
    script = open(scriptName, 'w')
    util.writeHeader(script, config, "sortbedgraph")
    script.write("sort -k1,1 -k2,2n " + "\\\n")
    script.write(inputDirectory + "/" + file + ".bedgraph > " + "\\\n")
    script.write(outputDirectory + "/" + file + ".bedgraph")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
