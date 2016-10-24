#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates stringtie scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="stringtie")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files.", default="../results/star2pass")
parser.add_argument("-o", "--outputDirectory", help="Output directory with stringtie results.", default="../results/stringtie")
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

trim = config.getboolean("project", "trim")
stranded = config.getboolean("project", "stranded")
genome = config.get("project", "genome")
genomeFile = config.get(genome, "genomeFile")
gtfFile = config.get(genome, "gtfFile")
processors = config.get("stringtie", "processors")

# Read samples file.
samplesFile = util.readSamplesFile()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

########################
# stringtie.sh scripts #
########################
for index, row in samplesFile.iterrows():
    sample = row["sample"]
    if not os.path.exists(os.path.join(outputDirectory, sample)):
        os.makedirs(os.path.join(outputDirectory, sample))
    scriptName = "stringtie_" + sample + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "stringtie")
    script.write("stringtie \\\n")
    script.write(os.path.relpath(os.path.join(inputDirectory, sample, sample + "Aligned.sortedByCoord.out.bam")) + " \\\n") 
    script.write("-G " + gtfFile + " \\\n")
    script.write("-o " + os.path.relpath(os.path.join(outputDirectory, sample, sample + ".gtf")) + " \\\n")
    script.write("-A " + os.path.relpath(os.path.join(outputDirectory, sample, sample + "_gene_abund.tab")) + " \\\n") 
    script.write("-C " + os.path.relpath(os.path.join(outputDirectory, sample, sample + "_cov_refs.gtf")) + " \\\n") 
    script.write("-p " + processors + " \\\n")
    script.write("-B" + " \\\n")
    script.write("&> " + scriptName + ".log")
if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
