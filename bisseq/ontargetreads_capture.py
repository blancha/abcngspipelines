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
parser = argparse.ArgumentParser(description="Generates scripts to count on primary target reads.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=ontargetreads_capture", default="ontargetreads_capture")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/bamtools_filter", default="../results/bamtools_filter")
parser.add_argument("-o", "--outputDirectory", help="Output directory with filterd BAM files. DEFAULT=../results/ontargetreads_capture", default="../results/ontargetreads_capture")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the scripts directory, cd to the scripts directory
util.cdMainScriptsDirectory()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
scriptsDirectory = os.path.abspath(args.scriptsDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

samples = util.getMergedsamples()

# Create script and output directories, if they do not exist yet.
util.makeDirectory(outputDirectory)
util.makeDirectory(scriptsDirectory)

# CD to scripts directory
os.chdir(scriptsDirectory)

# Write scripts
for sample in samples:
    scriptName =  "ontargetreads_capture_" + sample + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "ontargetreads_capture")
    # Reorder 
    
    script.write("bedtools intersect -bed ") + " \\\n")
    script.write("-abam" + os.path.join(inputDirectory, sample + ".filtered.bam") + " \\\n")
    script.write("-b " + "/sb/project/afb-431/BIF/SZY/SZY-BIF-P3/data/140501_RN4_PTSD3_RM_EPI_combined_coverage/140501_RN4_PTSD3_RM_EPI_capture_targets.bed"  + " \\\n")
    script.write("1> " + os.path.joint(outputDirectory, "intersect_" + sample + "_capture_targets.txt")
    script.write("2> " + scriptName + ".log")

    script.close()
