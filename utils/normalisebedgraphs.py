#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 21/03/2014

import argparse
import configparser
import glob
import os
import os.path
import pandas
import subprocess
import sys
import util

parser = argparse.ArgumentParser(description='Generate scripts to normalise bedgraph files by scaling factor.')
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=normalisebedgraph", default="normalisebedgraph")
parser.add_argument("-i", "--inputDirectory", help="Input directory with bedgraph files. DEFAULT=../results/bedgraph/ucsc", default="../results/bedgraph/ucsc")
parser.add_argument("-o", "--outputDirectory", help="Output directory with normalised bedgraph files. DEFAULT=../results/bedgraph/with_normalisation", default="../results/bedgraph/with_normalisation")
parser.add_argument("-f", "--scalingFactors", help="Tab delimited text file with scaling factors. DEFAULT=scaling_factors.txt", default="scaling_factors.txt")
parser.add_argument("-a", "--sampleSubDirectories", help="BAM files located in sample subdirectories in inputDirectory. DEFAULT=yes", choices=["yes", "no", "y", "n"], default="yes")
parser.add_argument("-e", "--extension", help="Extension identifying bedgraph files. DEFAULT=.bedgraph", default=".bedgraph") 
parser.add_argument("-t", "--stranded", help="Stranded bedgraph files. DEFAULT=yes", choices=["yes", "no", "y", "n"], default="yes")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
scalingFactors = args.scalingFactors
if (args.stranded == "yes") | (args.stranded == "y"):
    stranded = True
else:
    stranded = False

# Read scaling_factors file. Call script to create file first if it doesn't exist.
if not os.path.exists(args.scalingFactors):
    print("Scaling factors file does not exist. Calling scalingFactors.R")
    subprocess.call("scalingFactors.R", shell=True)
    # Exit program if file was not create.
    if not os.path.exists("scaling_factors.txt"):
        print("Error!\nscaling_factors.txt file could not be created.")
        print("Exiting program.")
        exit(1)

scalingFactors = pandas.read_csv(args.scalingFactors, sep="\t")

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(inputDirectory)

# Create script and output directories, if they do not exist yet.
util.makeDirectory(outputDirectory)
util.makeDirectory(scriptsDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
genome = config.get("project", "genome")
institute = config.get(genome, "institute")

# cd to scripts directory
os.chdir(scriptsDirectory)

for index, row in scalingFactors.iterrows():
    sample = row["name"]
    scalingFactor = str(row["scaling_factor"])
    # Create script file.
    scriptName = 'normalisebedgraph_' + sample + '.sh'
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "normalisebedgraph_batch")
    script.write("normalisebedgraph.py " + "\\\n")
    script.write("--before_normalisation " + os.path.relpath(os.path.join(inputDirectory, sample + ".bedgraph")) + " \\\n")
    script.write("--after_normalisation " + os.path.relpath(os.path.join(outputDirectory, sample + ".bedgraph")) + " \\\n")
    script.write("--scaling_factor " + scalingFactor + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.close()
    if stranded: 
        # Create positive strand script file.
        scriptName = 'normalisebedgraph_' + sample + '_positive.sh'
        script = open(scriptName, 'w')
        if header:
            util.writeHeader(script, config, "normalisebedgraph_batch")
        script.write("normalisebedgraph.py " + "\\\n")
        script.write("--before_normalisation " + os.path.relpath(os.path.join(inputDirectory, sample + "_positive.bedgraph")) + " \\\n")
        script.write("--after_normalisation " + os.path.relpath(os.path.join(outputDirectory, sample + "_positive.bedgraph")) + " \\\n")
        script.write("--scaling_factor " + scalingFactor + " \\\n")
        script.write("&> " + scriptName + ".log")
        script.close()
        # Create negative strand script file.
        scriptName = 'normalisebedgraph_' + sample + '_negative.sh'
        script = open(scriptName, 'w')
        if header:
            util.writeHeader(script, config, "normalisebedgraph_batch")
        script.write("normalisebedgraph.py " + "\\\n")
        script.write("--before_normalisation " + os.path.relpath(os.path.join(inputDirectory, sample + "_negative.bedgraph")) + " \\\n")
        script.write("--after_normalisation " + os.path.relpath(os.path.join(outputDirectory, sample + "_negative.bedgraph")) + " \\\n")
        script.write("--scaling_factor " + scalingFactor + " \\\n")
        script.write("&> " + scriptName + ".log")
        script.close()
if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
