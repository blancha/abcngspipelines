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
parser = argparse.ArgumentParser(description="Generates bamtools split scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="splitbam")
parser.add_argument("-e", "--extension", help="Extension to identify BAM files to sort. DEFAULT=.bam", default=".bam")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/bsmap", default="../results/bsmap")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the scripts directory, cd to the scripts directory
util.cdMainScriptsDirectory()

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

# Get samples and conditions
samples = util.getMergedsamples()
samples = sorted(samples, reverse=True) # To always put wt first, put the list in reverse alpabetical order

# Create scripts directory, if it does not exist yet, and cd to it.
util.makeDirectory(scriptsDirectory)
os.chdir(scriptsDirectory)

for sample in samples:
    # Create script
    scriptName = "splitbam_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "samtoolsIndex")
    script.write("bamtools split " + "\\\n")    
    script.write("-tag ZS " + "\\\n")
    inputFile = glob.glob(inputDirectory + "/" + sample + "*" + args.extension)[0]
    script.write("-in " + inputFile + " \\\n")
    script.write("&> " + scriptName + ".log")

    script.write("\n\n")

    inputFile_basename = os.path.join(inputDirectory, os.path.splitext(os.path.basename(inputFile))[0])

    script.write("bamtools merge " + "\\\n")
    script.write("-in " +  inputFile_basename + ".TAG_ZS_++.bam " + "\\\n")
    script.write("-in " +  inputFile_basename + ".TAG_ZS_-+.bam " + "\\\n")
    script.write("-out " + inputFile_basename + ".top.bam " + "\\\n")
    script.write("&>> " + scriptName + ".log")

    script.write("\n\n")

    script.write("bamtools merge " + "\\\n")
    script.write("-in " + inputFile_basename + ".TAG_ZS_-+.bam " + "\\\n")
    script.write("-in " + inputFile_basename + ".TAG_ZS_--.bam " + "\\\n")
    script.write("-out " + inputFile_basename + ".bottom.bam "  + "\\\n")
    script.write("&>> " + scriptName + ".log")

    script.write("\n\n")

    script.write("samtools sort " + "\\\n")
    script.write(inputFile_basename + ".top.bam " + "\\\n")
    script.write(inputFile_basename + ".top.sorted " + "\\\n")
    script.write("&>> " + scriptName + ".log")

    script.write("\n\n")

    script.write("samtools sort " + "\\\n")
    script.write(inputFile_basename + ".bottom.bam " + "\\\n")
    script.write(inputFile_basename + ".bottom.sorted " + "\\\n")
    script.write("&>> " + scriptName + ".log")

script.close()

