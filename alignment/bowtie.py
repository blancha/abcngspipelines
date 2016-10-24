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
parser = argparse.ArgumentParser(description="Generates Bowtie2 scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="bowtie")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files.", default="../data/FASTQ_files/untrimmed/")
parser.add_argument("-o", "--outputDirectory", help="Output directory with bowtie results.", default="../results/bowtie/")
#parser.add_argument("-x", "--index", help="Index BAM files.", choices=["yes", "no", "y", "n"], default="yes")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
#if (args.index.lower() == "yes") | (args.index.lower() == "y"):
#    index=True
#else:
#    index=False

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
genome = config.get("project", "genome")
gtfFile = config.get(genome, "gtfFile")
bowtieIndex = config.get(genome, "bowtie2Index")
processors = config.get('bowtie', 'processors')
unaligned = config.getboolean("bowtie", "un")
no_discordant = config.getboolean("bowtie", "no_discordant")

# Read samples file.
samplesFile = util.readSamplesFile(inputDirectory)

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directories, if they do not exist yet..
if not os.path.exists(outputDirectory): 
    os.makedirs(outputDirectory)

# Cycle through all the samples and write the bowtie scripts.
for index, row in samplesFile.iterrows():
    sample = row["sample"]
    # Create output directories
    if not os.path.exists(outputDirectory + "/" + sample):
        os.mkdir(outputDirectory +"/" + sample)
    file_R1 = row["file_r1"]
    if "file_r2" in samplesFile.columns:
        file_R2 = row["file_r2"]
    # Create script file.
    scriptName = 'bowtie_' + sample + '.sh'
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "bowtie")
    script.write("bowtie2 " + "\\\n")
    script.write("-p " + processors + " \\\n");
    if no_discordant:
       script.write("--no-discordant" + " \\\n")
    script.write("-x " + bowtieIndex + " \\\n")
    if "file_r2" in samplesFile.columns:
        script.write("-1 " + os.path.relpath(os.path.join(inputDirectory, file_R1)) + " \\\n")
        script.write("-2 " + os.path.relpath(os.path.join(inputDirectory, file_R2)) + " \\\n")
    else:
        script.write("-U " + os.path.relpath(os.path.join(inputDirectory, file_R1)) + " \\\n")
    if unaligned:
        script.write("--un " + os.path.relpath(os.path.join(outputDirectory, sample)) + " \\\n")
    script.write("1> " +  os.path.relpath(os.path.join(outputDirectory, sample, sample + ".sam")) + " \\\n")
    script.write("2> " + scriptName + ".log") 
    script.write("\n\n")

    script.write("samtools view -hbS " + os.path.relpath(os.path.join(outputDirectory, sample, sample + ".sam ")) + "\\\n")
    script.write("1> " + os.path.relpath(os.path.join(outputDirectory, sample, sample + ".bam")) + " \\\n")
    script.write("2>> " + scriptName + ".log") 
    script.write("\n\n")

    script.write("rm " + os.path.relpath(os.path.join(outputDirectory, sample, sample + ".sam")) + " \\\n")
    script.write("2>> " + scriptName + ".log")    

    #if index: 
    #    script.write("\n\n")
    #    script.write("samtools index " + os.path.relpath(os.path.join(outputDirectory, sample, sample + ".bam")) + " \\\n")
    #    script.write("2>> " + scriptName + ".log")    

    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
