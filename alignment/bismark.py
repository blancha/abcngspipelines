#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates Bismark scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=bismark", default="bismark")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files. DEFAULT=../data/FASTQ_files/trimmed", default="../data/FASTQ_files/trimmed")
parser.add_argument("-o", "--outputDirectory", help="Output directory with Tophat results. DEFAULT=../results/bismark", default="../results/bismark")
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

genome = config.get("project", "genome")
bisulfiteGenomeFolder = config.get(genome, "bisulfiteGenomeFolder")

# Read samples file.
samplesFile = util.readsamplesFile()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Cycle through all the samples and write the Bismark scripts.
for index, row in samplesFile.iterrows():
    sample = row["sample"]
    if "Lane" in samplesFile.columns:
        sample= sample + "_lane_" + str(row["Lane"])
    # Create output directories
    if not os.path.exists(outputDirectory + "/" +  sample):
        os.mkdir(outputDirectory + "/" + sample)    
    file_R1 = row["File_R1"]
    # If trimmed with trim_galore, extensions have changed.
    file_R1 = file_R1.replace("R1.fastq.gz", "R1_val_1.fq.gz")
    if "File_R2" in samplesFile.columns:
        file_R2 = row["File_R2"]
        # If trimmed with trim_galore, extensions have changed.
        file_R2 = file_R2.replace("R2.fastq.gz", "R2_val_2.fq.gz")
    # Create script file.
    scriptName = 'bismark_' + sample + '.sh'
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "bismark")
    script.write("bismark" + " \\\n")
    script.write("--bowtie2" + " \\\n")
    script.write("--basename " + sample + " \\\n")
    script.write("-output_dir " + outputDirectory + "/" + sample + " \\\n")
    script.write(bisulfiteGenomeFolder + " \\\n")
    script.write("-1 " + inputDirectory + "/" + file_R1 + " \\\n")
    script.write("-2 " + inputDirectory + "/" + file_R2 + " \\\n")
    script.write("&> " + scriptName + ".log")    

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
