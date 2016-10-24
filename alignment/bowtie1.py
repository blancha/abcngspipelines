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
parser = argparse.ArgumentParser(description="Generates Bowtie1 scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="bowtie")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files.", default="../data/FASTQ_files/trimmed/")
parser.add_argument("-o", "--outputDirectory", help="Output directory with bowtie results.", default="../results/bowtie/")
parser.add_argument("-g", "--genome", help="Reference genome. If left unspecified, reference genome will be read from configuration file DEFAULT=Configuration file value", default="Configuration file value")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
genome = args.genome

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
stranded = config.getboolean("project", "stranded")
if genome is "Configuration file value":
    if config.has_option("project", "genome"):
        genome = config.get("project", "genome")
    else:
        sys.exit("No reference genome was specified at the command line or in the configuration file.") 
genomeFile = config.get(genome, "genomeFile")
gtfFile = config.get(genome, "gtfFile")
bowtieIndex = config.get(genome, "bowtie1Index")
processors = config.get("bowtie1", 'processors')
k = config.get("bowtie1", "k")
v = config.get("bowtie1", "v")

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
    file_r1 = row["file_r1"]
    if "file_r2" in samplesFile.columns:
        file_R2 = row["file_r2"]
    # Create script file.
    scriptName = 'bowtie_' + sample + '.sh'
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "bowtie")
    script.write("bowtie " + "\\\n")
    script.write("-S" + " \\\n")
    script.write("-p " + processors + " \\\n")
    script.write("-k " + k + " \\\n")
    script.write("-v " + v + " \\\n")
    script.write(bowtieIndex + " \\\n")
    script.write(os.path.relpath(os.path.join(inputDirectory, file_r1)) + " \\\n")
    script.write("1> " + os.path.relpath(os.path.join(outputDirectory, sample, sample + ".sam")) + " \\\n")
    script.write("2> " + scriptName + ".log") 
#    script.write("\n\n")

#    script.write("samtools view -hbS " + outputDirectory + "/" + sample + "/" + sample + ".sam " + "\\\n")
#    script.write("1> " + outputDirectory + "/" + sample + "/" + sample + ".bam " + "\\\n")
#    script.write("2>> " + scriptName + ".log") 
#    script.write("\n\n")

#    script.write("rm " + outputDirectory + "/" + sample + "/" + sample + ".sam" + " \\\n")
#    script.write("2>> " + scriptName + ".log")    

    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
