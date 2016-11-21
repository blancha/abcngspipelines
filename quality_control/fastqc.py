#!/usr/bin/env python3

# Version 1.2
# Author Alexis Blanchet-Cohen
# Last modified: 27/10/2016

import argparse
import glob
import os
import os.path
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description='Generates FastQC scripts.')
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=fastqc", default="fastqc")
parser.add_argument("-o", "--outputDirectory", help="Output directory with FastQC results. DEFAULT=../results/fastqc", default="../results/fastqc")
parser.add_argument("-l", "--laneDirectories", help="Create a different directory for each lane. Useful when the FASTQ files for different lanes have the same name.", choices=["yes", "no", "y", "n"], default="no")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Read samples file.
samplesDataFrame = util.readSamplesFile()
samples = samplesDataFrame["sample"].tolist()

# Create script and output directories, if they do not exist yet.
util.makeDirectory(outputDirectory, recursive=True)
util.makeDirectory(scriptsDirectory, recursive=True)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
use_modules = config.getboolean("server", "use_modules")
if use_modules:
    module = config.get("fastqc", "module")

# Change to scripts directory
os.chdir(scriptsDirectory)

# Cycle through all the samples and write the fastqc scripts.
for index, row in samplesDataFrame.iterrows():
    sample = row["sample"]
    file_r1 = row["file_r1"]
    file_r2 = row["file_r2"]
    files = [file_r1, file_r2]
    if "lane" in samplesDataFrame.columns:
        lane = "_" + row["lane"]
    else:
        lane = ""
    # Write script
    for file in files:                  
        scriptName = "fastqc_" + os.path.basename(file) +lane + ".sh"
        script = open(scriptName, "w")
        if header:
            util.writeHeader(script, config, "fastqc")
            if use_modules:
                util.write("module load " + module + "\\\n\\\n")    
            script.write("fastqc " + "\\\n")
            if (args.laneDirectories == "no") | (args.laneDirectories == "no"):
                script.write("--outdir " + os.path.relpath(outputDirectory)  + " \\\n")
            elif "lane" in samplesDataFrame.columns:
                if not os.path.exists(os.path.join(outputDirectory, row["lane"])):
                    os.mkdir(outputDirectory, row["lane"])
                script.write("--outdir " + os.path.relpath(outputDirectory, row["lane"])  + " \\\n")
            script.write(file + " \\\n")
            script.write("&> " + scriptName + ".log")
            script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)

