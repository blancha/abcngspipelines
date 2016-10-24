#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 04/03/2014

import argparse
import configparser
import glob
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Concatenates all SRA FASTQ files from the same sample.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="catfastqfilessra")
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

samplesFile = pandas.read_csv(os.path.join(inputDirectory, "SraRunTable.txt"), sep="\t")

# Keep only the unique sample names
unique_samples = samplesFile.drop_duplicates(subset="BioSample_s").BioSample_s.tolist()

# Write the script
scriptName = 'catfastqfilessra.sh'
script = open(scriptName, 'w')    
if header:
    util.writeHeader(script, config, "catfastqfilessra")
    first = True
    for unique_sample in unique_samples:
        runs = samplesFile[samplesFile["BioSample_s"] == unique_sample].Run_s.tolist()   
        for end in ["1", "2"]:
            script.write("cat " + "\\\n")
            for run in runs:
                script.write(os.path.relpath(os.path.join(inputDirectory, run)) + "_" + end + ".fastq" + " \\\n")
            script.write("1> " + os.path.relpath(os.path.join(outputDirectory, unique_sample)) + "_" + end + ".fastq" + " \\\n")
            if(first==True):
                script.write("2> " + scriptName + ".log")
                first=False
            else:
                script.write("2>> " + scriptName + ".log")
            script.write("\n\n")
script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
