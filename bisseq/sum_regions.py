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
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=sum_regions", default="sum_regions")
parser.add_argument("-i", "--inputFile", help="Target BED file. DEFAULT=../results/target_intervals", default="../../data/140501_RN4_PTSD3_RM_EPI_combined_coverage/140501_RN4_PTSD3_RM_EPI_capture_targets.bed")
parser.add_argument("-o", "--outputDirectory", help="Output directory with padded target BED file. DEFAULT=../results/padded_target", default="padded_target")
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
scriptName =  "sum_regions.sh"
script = open(scriptName, "w")
if header:
    util.writeHeader(script, config, "sum_regions")
    
script.write("bedtools genomecov " + "\\\n")
script.write("-i " + inputFile + " \\\n")
script.write("-g " + chromInfo + " -max 1 | " + " \\\n")
script.write("grep -P \"genome\\t1\" | cut -f 3" + "\\\n")
script.write("&> " + scriptName + ".log")

script.close()
