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
parser = argparse.ArgumentParser(description="Generates Picard target list scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=add_padding_targets", default="add_padding_targets")
parser.add_argument("-i", "--inputFile", help="Target BED file. DEFAULT=../results/target_intervals", default="../../data/140501_RN4_PTSD3_RM_EPI_combined_coverage/140501_RN4_PTSD3_RM_EPI_capture_targets.bed")
parser.add_argument("-o", "--outputDirectory", help="Output directory with padded target BED file. DEFAULT=../results/padded_targets", default="../results/padded_targets")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the scripts directory, cd to the scripts directory
util.cdMainScriptsDirectory()

# Process the command line arguments.
outputDirectory = os.path.abspath(args.outputDirectory)
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
inputFile = args.inputFile

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

genome = config.get("project", "genome")
chromInfo = config.get(genome, "chromInfo")

# Create script and output directories, if they do not exist yet.
util.makeDirectory(outputDirectory)
util.makeDirectory(scriptsDirectory)

# CD to scripts directory
os.chdir(scriptsDirectory)

# Write scripts
scriptName =  "add_padding_targets.sh"
script = open(scriptName, "w")
if header:
    util.writeHeader(script, config, "add_padding_targets")
    
script.write("bedtools slop " + "\\\n")
script.write("-i " + inputFile + " \\\n")
script.write("-b 100 -g " + chromInfo + " | " + " \\\n")
script.write("bedtools sort -i - | " + "\\\n")
script.write("bedtools merge -i - " + "\\\n")
script.write("1> " + os.path.join(outputDirectory, "140501_RN4_PTSD3_RM_EPI_padded_capture_target.bed") + " \\\n")
script.write("2> " + scriptName + ".log")

script.close()
