#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates Tophat scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=tophat", default="tophat")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files. DEFAULT=../data/FASTQ_files/untrimmed", default="../data/FASTQ_files/untrimmed")
parser.add_argument("-o", "--outputDirectory", help="Output directory with Tophat results. DEFAULT=../results/tophat", default="../results/tophat")
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

trim = config.getboolean("project", "trim")
stranded = config.getboolean("project", "stranded")
genome = config.get("project", "genome")
genomeFile = config.get(genome, "genomeFile")
gtfFile = config.get(genome, "gtfFile")
bowtie2Index = config.get(genome, "bowtie2Index")
processors = config.get('tophat', 'processors')
no_novel_juncs = config.getboolean("tophat", "no_novel_juncs")
no_discordant = config.getboolean("tophat", "no_discordant")

# Read samples file.
samplesFile = util.readSamplesFile()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Cycle through all the samples and write the tophat scripts.
for index, row in samplesFile.iterrows():
    sample = row["sample"]
    if "lane" in samplesFile.columns:
        sample= sample + "_lane_" + str(row["lane"])
    # Create output directories
    if not os.path.exists(outputDirectory + "/" +  sample):
        os.mkdir(outputDirectory + "/" + sample)    
    file_R1 = row["file_r1"]
    if "file_r2" in samplesFile.columns:
        file_R2 = row["file_r2"]
    # Create script file.
    scriptName = "tophat_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "tophat")
    script.write("#Deactivate Python 3 virtual environment, and activate Python 2 virtual environment to be able to use TopHat.\n")
    virtual_env_directory = os.path.dirname(os.environ["VIRTUAL_ENV"])
    script.write("source " + os.path.join(virtual_env_directory, "python2.7/bin/activate") + "\n\n")    
    script.write("tophat " + "\\\n")
    script.write("--rg-library \"L\" --rg-platform \"ILLUMINA\" " + "\\\n")
    script.write("--rg-platform-unit \"X\" --rg-sample \"" + sample + "\" --rg-id \"runX\" " + "\\\n")
    if no_novel_juncs:
        script.write("--no-novel-juncs" + " \\\n")
    if no_discordant:
        script.write("--no-discordant" + " \\\n")
    script.write("-p " + processors + " \\\n"); 
    if stranded:
        script.write("--library-type fr-firststrand " + "\\\n")
    script.write("-G " + gtfFile + " \\\n") 
    script.write("-o " + os.path.relpath(os.path.join(outputDirectory, sample)) + " \\\n")
    script.write(bowtie2Index + " \\\n")
    if trim:
        if "untrimmed" in inputDirectory:
            print("Since FASTQ files are trimmed, replacing untrimmed input directory with trimmed input directory")
            inputDirectory = inputDirectory.replace("untrimmed", "trimmed")  
    script.write(inputDirectory + "/" + file_R1 + " \\\n")
    if "file_r2" in samplesFile.columns:
        script.write(inputDirectory + "/" + file_R2 + " \\\n")
    script.write("&> " + scriptName + ".log")    

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
