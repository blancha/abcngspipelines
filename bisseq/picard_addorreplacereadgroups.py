#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import subprocess
import sys
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates picard_addorreplacereadgroups scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="picard_addorreplacereadgroups")
parser.add_argument("-e", "--extension", help="Extension to identify SAM files which must be converted to BAM files. DEFAULT=.sam", default=".sam")
parser.add_argument("-i", "--inputDirectory", help="Input directory with SAM files. DEFAULT=../results/bsmap", default="../results/bsmap")
parser.add_argument("-o", "--outputDirectory", help="Output directory with sorted SAM files. DEFAULT=..results/bsmap", default="../results/bsmap")
#parser.add_argument("-u", "--subDirectories", help="samples in subdirectories. DEFAULT=..yes", choices=["yes", "no", "y", "n"], default="yes")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the scripts directory, cd to the scripts directory
util.cdMainScriptsDirectory()

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

picard_folder = config.get("picard", "folder")

# Get samples and conditions
samples = util.getMergedsamples()

# Create scripts directory, if it does not exist yet, and cd to it.
util.makeDirectory(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
util.makeDirectory(outputDirectory)

for sample in samples:
    # Create script
    scriptName = "picard_addorreplacereadgroups_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "picard_addorreplacereadgroups")
    script.write("java -jar " + picard_folder + "AddOrReplaceReadGroups.jar " + "\\\n")
    script.write("VALIDATION_STRINGENCY=LENIENT " + "\\\n")
    inputFile = glob.glob(inputDirectory + "/" + sample + args.extension)[0]
    script.write("INPUT=" + inputFile + " \\\n")
    outputFile = outputDirectory + "/" + sample + ".bam"
    script.write("OUTPUT=" + outputFile + " \\\n")
    script.write("RGID=" + sample + " \\\n")
    script.write("RGLB=" + sample + " \\\n")
    script.write("RGPL=Illumina" + " \\\n")
    script.write("RGSM=" + sample  + " \\\n")
    script.write("RGPU=platform_unit" + " \\\n")
    script.write("&> " + scriptName + ".log")
script.close()

