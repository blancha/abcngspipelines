#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import os.path
import subprocess
import sys
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates samtools index scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Directory in which scripts will be placed. Directory will be created, if it doesn't exist yet. DEFAULT=samtools_index", default="samtools_index")
parser.add_argument("-e", "--extension", help="Extension to identify BAM files to index. DEFAULT=.bam", default=".bam")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=..results/tophat", default="../results/tophat")
parser.add_argument("-a", "--sampleSubDirectories", help="BAM files located in sample subdirectories in inputDirectory. DEFAULT=yes", choices=["yes", "no", "y", "n"], default="yes")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# If not in the scripts directory, cd to the scripts directory
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
extension = args.extension

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

processors = config.get("samtools", "processors")

# Read samples file.
samplesFile = util.readSamplesFile()
samples = samplesFile["sample"].tolist()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# List BAM files
if (args.sampleSubDirectories == "yes") | (args.sampleSubDirectories == "y"):
    files = []
    samples = os.listdir(inputDirectory)
    for sample in samples:
        cwd = os.getcwd()
        os.chdir(os.path.join(inputDirectory, sample))
        sample_files = glob.glob("*" + extension)
        files.append(sample_files[0])
        os.chdir(cwd)
else:
    samples = files = os.listdir(inputDirectory)

for sample,file in zip(samples, files):
    # Create script
    scriptName = "samtools_index_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "samtools_index")
    script.write("samtools index " + "\\\n")
    inputFile = os.path.join(inputDirectory, sample, file)
    script.write(inputFile + " \\\n")
    script.write("&> " + scriptName + ".log")

script.close()

