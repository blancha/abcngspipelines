#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Renames accepted_hits.bam to sample_name.bam scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=tophat", default="rename_and_index_accepted_hits_bam")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files. DEFAULT=../data/FASTQ_files/untrimmed", default="../results/tophat")
parser.add_argument("-o", "--outputDirectory", help="Output directory with renamed BAM files DEFAULT=../results/tophat", default="../results/tophat")
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

# Read samples file. Create file first if it doesn't exist.
samplesFile = pandas.read_csv("samples.txt", sep="\t")

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Cycle through all the samples and write the tophat scripts.
for index, row in samplesFile.iterrows():
    sample = row["Sample"]
    if "Lane" in samplesFile.columns:
        sample= sample + "_lane_" + str(row["Lane"])
    # Create output directories
    if not os.path.exists(outputDirectory + "/" +  sample):
        os.mkdir(outputDirectory + "/" + sample)    
    scriptName = "rename_and_index_accepted_hits_bam_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "rename_and_index_accepted_hits_bam")
    script.write("mv " + "\\\n")
    script.write(os.path.relpath(os.path.join(inputDirectory, sample, "accepted_hits.bam")) + " \\\n")
    script.write(os.path.relpath(os.path.join(outputDirectory, sample, sample + ".bam")) + " \\\n")
    script.write("&> " + scriptName + ".log")    
    script.write("\n\n")
    script.write("samtools index " + os.path.relpath(os.path.join(outputDirectory, sample, sample + ".bam")) + " \\\n")
    script.write("&>> " + scriptName + ".log")    

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
