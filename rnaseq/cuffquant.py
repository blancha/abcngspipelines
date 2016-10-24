#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
from collections import OrderedDict
import glob
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates Cuffquant scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="cuffquant")
parser.add_argument("-e", "--extension", help="Extension to identify BAM files to sort. DEFAULT=.bam", default=".bam")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files.", default="../results/star")
parser.add_argument("-o", "--outputDirectory", help="Output directory with Cuffquant results.", default="../results/cuffquant")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
extension = args.extension

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

processors = config.get("cuffquant", "processors")
stranded = config.getboolean("project", "stranded")
genome = config.get("project", "genome")
genomeFile = config.get(genome, "genomeFile")
gtfFile = config.get(genome, "gtfFile")
if config.has_option("cuffquant", "maskFile"):
    maskFile = config.get("cuffquant", "maskFile")

# Read samples file. Create file first if it doesn't exist.
if not os.path.exists("samples.txt"):
    subprocess.call("samples.py", shell=True)
samplesFile = pandas.read_csv("samples.txt", delim_whitespace=True)

# Get samples
samples = samplesFile["sample"]

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directories, if they have not been created yet.
for sample in samples:
    if not os.path.exists(os.path.join(outputDirectory, sample)):
        os.makedirs(os.path.join(outputDirectory, sample))

# Write the cuffquant scripts.
for sample in samples:
    scriptName = "cuffquant_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "cuffquant")
    script.write("cuffquant" + " \\\n")
    script.write("-p " + processors + " \\\n")
    script.write("--no-effective-length-correction " + "\\\n")
    if stranded:
        script.write("--library-type fr-firststrand" + " \\\n")
    if config.has_option("cuffquant", "maskFile"):
        script.write("--mask-file " + os.path.relpath(maskFile) + " \\\n")
    script.write("-u -b " + genomeFile + " \\\n") 
    script.write("-o " + os.path.relpath(os.path.join(outputDirectory, sample)) + " \\\n")
    script.write(gtfFile + " \\\n")
    bamFile = glob.glob(os.path.join(inputDirectory, sample, "*" + extension))[0] 
    script.write(os.path.relpath(os.path.join(inputDirectory, sample, bamFile)) + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
