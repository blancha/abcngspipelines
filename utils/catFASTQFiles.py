#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 04/03/2014

import argparse
import configparser
import glob
import os
import os.path
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Concatenates all FASTQ files from the same lane..")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="catFASTQFiles")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files. DEFAULT=../data/FASTQ_files/untrimmed", default="../data/FASTQ_files/untrimmed")
parser.add_argument("-o", "--outputDirectory", help="Output directory with concatenated FASTQ files ../data/FASTQ_files/untrimmed.", default="../data/FASTQ_files/untrimmed")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.mkdir(outputDirectory)


############################
# Group the files by lane. #
############################
files = glob.glob(os.path.join(inputDirectory, "*fastq")) + glob.glob(os.path.join(inputDirectory, "*fastq.gz"))

# Identify all unique files, regardless of lane.
samples = []
for file in files:
    fields = file.split(".")
    samples.append(".".join(fields[-3:]))
# Remove duplicates, and sort lanes.    
samples = sorted(list(set(samples)))

all_samples_grouped_by_lane = []
for sample in samples:
    one_sample_grouped_by_lane = sorted(glob.glob(inputDirectory + "/*" + sample))
    all_samples_grouped_by_lane.append(one_sample_grouped_by_lane)

# Write the script
scriptName = 'catFASTQFiles.sh'
script = open(scriptName, 'w')    
if header:
    util.writeHeader(script, config, "catFASTQFiles")

for sample, one_sample_grouped_by_lane in zip(samples, all_samples_grouped_by_lane):
    script.write("cat " + "\\\n")
    script.write(" \\\n".join(one_sample_grouped_by_lane) + " \\\n")
    script.write("1> " + os.path.relpath(os.path.join(outputDirectory, sample)) + " \\\n")
    script.write("2>> " + scriptName + ".log")
    script.write("\n\n")

script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
