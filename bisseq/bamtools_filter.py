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
parser = argparse.ArgumentParser(description="Generates Picard tools MarkDuplicates scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=bamtools_filter", default="bamtools_filter")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/bamtools_merge", default="../results/bamtools_merge")
parser.add_argument("-o", "--outputDirectory", help="Output directory with filterd BAM files. DEFAULT=../results/bamtools_filter", default="../results/bamtools_filter")
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
    
    scriptName =  "bamtools_filter_" + sample + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "bamtools_filter")
    # Reorder 
    script.write("bamtools filter " + "\\\n")
    script.write("-isMapped true" + " \\\n")
    script.write("-isPaired true" + " \\\n")
    script.write("-isProperPair true" + " \\\n")
    script.write("-forceCompression" + " \\\n")
    script.write("-in " + os.path.join(inputDirectory, sample + ".rmdups.bam") + " \\\n")
    script.write("-out " + os.path.join(outputDirectory, sample + ".filtered.bam") + " \\\n")
    script.write("&> " + scriptName + ".log")

    script.close()
