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
parser = argparse.ArgumentParser(description="Generates htseqcount scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="htseqcount")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files.", default="../results/tophat")
parser.add_argument("-o", "--outputDirectory", help="Output directory with htseq-count results.", default="../results/htseqcount")
parser.add_argument("-e", "--extension", help="Extension of BAM files on which to run htseq-count", default=".bam")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
extension = args.extension

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

trim = config.getboolean("project", "trim")
stranded = config.getboolean("project", "stranded")
genome = config.get("project", "genome")
genomeFile = config.get(genome, "genomeFile")
gtfFile = config.get(genome, "gtfFile")

# Read samples file.
samplesFile = util.readSamplesFile()
samples = samplesFile["sample"]

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

#########################
# htseqcount.sh scripts #
#########################
for sample in samples:
    scriptName = "htseqcount_" + sample + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "htseqcount")
    script.write("#Deactivate Python 3 virtual environment, and activate Python 2 virtual environment to be able to use TopHat.\n")
    virtual_env_directory = os.path.dirname(os.environ["VIRTUAL_ENV"])
    script.write("source " + os.path.join(virtual_env_directory, "python2.7/bin/activate") + "\n\n") 
    script.write("samtools sort -n" + " \\\n")
    script.write("-o " + os.path.relpath(os.path.join(inputDirectory, sample, sample + "_sorted_by_read_name" + extension)) + " \\\n")
    script.write(os.path.relpath(os.path.join(inputDirectory, sample, sample + extension)) + " \\\n")
    script.write("&> " + scriptName + ".log")

    script.write("\n\n")

    script.write("htseq-count \\\n")
    script.write("--format=bam" + " \\\n")
    if stranded:
        script.write("--stranded=reverse" + " \\\n")
    else:
        script.write("--stranded=no" + " \\\n")
    script.write("--mode=intersection-strict" + " \\\n")
    script.write(os.path.relpath(os.path.join(inputDirectory, sample, sample + "_sorted_by_read_name" + extension)) + " \\\n")
    script.write(gtfFile + " \\\n")
    script.write("1> " + outputDirectory + "/" + sample + ".txt " + "\\\n")
    script.write("2>> " + scriptName + ".log") 

    script.write("\n\n")

    script.write("rm " + os.path.relpath(os.path.join(inputDirectory, sample, sample + "_sorted_by_read_name" + extension)) + " \\\n")
    script.write("&>> " + scriptName + ".log")

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
