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
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=target_interval", default="target_interval")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/bamtools_merge", default="../results/bamtools_filter")
parser.add_argument("-o", "--outputDirectory", help="Output directory with filterd BAM files. DEFAULT=../results/target_interval", default="../results/target_interval")
parser.add_argument("-r", "--target_BED", help="Target BED file. DEFAULT=/sb/project/afb-431/BIF/SZY/SZY-BIF-P3/data/140501_RN4_PTSD3_RM_EPI_combined_coverage/140501_RN4_PTSD3_RM_EPI_primary_targets.bed", default="/sb/project/afb-431/BIF/SZY/SZY-BIF-P3/data/140501_RN4_PTSD3_RM_EPI_combined_coverage/140501_RN4_PTSD3_RM_EPI_primary_targets.bed")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the scripts directory, cd to the scripts directory
util.cdMainScriptsDirectory()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
target_BED = args.target_BED

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
    scriptName =  "target_interval_" + sample + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "target_interval")
    # Reorder 
    
    script.write("samtools view -H" + " \\\n")
    script.write(os.path.join(inputDirectory, sample + ".filtered.bam") + " \\\n")
    script.write("1> " + os.path.join(outputDirectory, sample + "_bam_header.txt") + " \\\n")
    script.write("2> " + scriptName + ".log")

    script.write("\n\n")

    script.write("cat " + target_BED + " | " + "\\\n")
    script.write("gawk '{print $1 \"\\t\" $2+1 \"\\t\" $3 \"\\t+\\tinterval_\" NR}' " + "\\\n")
    script.write("1> " + os.path.join(outputDirectory, sample + "_design_target_body.txt") + " \\\n")
    script.write("2>> " + scriptName + ".log")

    script.write("\n\n")

    script.write("cat " + "\\\n")
    script.write(os.path.join(outputDirectory, sample + "_bam_header.txt") + " \\\n")
    script.write(os.path.join(outputDirectory, sample + "_design_target_body.txt") + " \\\n")
    script.write("1> " + os.path.join(outputDirectory, sample + "_design_target_intervals.txt") + " \\\n")
    script.write("2>> " + scriptName + ".log")

    script.close()
