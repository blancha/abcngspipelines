#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 12/04/2014

import argparse
import glob
import os
import os.path
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates HOMER makeTagDirectory scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=homer_makeTagDirectory", default="homer_makeTagDirectory")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files. DEFAULT=../results/bowtie", default="../results/bowtie")
parser.add_argument("-o", "--outputDirectory", help="Output directory with Tophat results. DEFAULT=../results/homer_makeTagDirectory", default="../results/homer_makeTagDirectory")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

genome = config.get('project', 'genome')

samples = util.getSamples()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Generate tag directories for all the BAM files.
for sample in samples:
    inputDirectory_sample = os.path.join(inputDirectory, sample)
    # Create output directories
    outputDirectory_sample = os.path.join(outputDirectory, sample)
    if not os.path.exists(outputDirectory_sample):
        os.makedirs(outputDirectory_sample)
    # Create script file
    scriptName = 'makeTagDirectory_' + sample + '.sh'
    script = open(scriptName, 'w')
    util.writeHeader(script, config, "makeTagDirectory")
    script.write("makeTagDirectory " + "\\\n")
    script.write(os.path.relpath(os.path.join(outputDirectory, sample)) + " \\\n")
    script.write(os.path.relpath(os.path.join(inputDirectory, sample, sample + "_sorted_by_coordinates.bam")) + " \\\n")
    script.write("-illuminaPE -checkGC -genome " + genome + " \\\n")
    script.write("&> " + scriptName + ".log")    
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
