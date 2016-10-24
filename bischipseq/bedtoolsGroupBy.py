#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 21/05/2014

import argparse
import glob
import os
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description='Bedtools groupby narrowPeaks with methylation calls.')
parser.add_argument("-c", "--scriptsDirectory", help="Scripts directory.", default="bedtoolsGroupBy")
parser.add_argument("-m", "--inputDirectory", help="Input directory with BED files.", default="../bedtoolsIntersect/narrowPeaks_and_cov/")
parser.add_argument("-o", "--outputDirectory", help="Output directory with grouped BED files.", default="../bedtoolsGroupBy/narrowPeaks_and_cov/")
parser.add_argument("-s", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

samples = util.getMergedsamples()

# Read configuration files.
config = util.readConfigurationFiles()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

for sample in samples:
    # Create script file.
    scriptName = 'bedtoolsGroupBy_' + sample + '.sh'
    script = open(scriptName, 'w')
    util.writeHeader(script, config, "bedtoolsGroupBy")
    script.write("bedtools groupby " + "\\\n")
    script.write("-i " + inputDirectory +  "/" + sample + ".bed " + "\\\n")
    script.write("-g 1,2,3,4 " + "\\\n")
    script.write("-o mean,sum,sum " + "\\\n")
    script.write("-opCols 14,15,16 " + "\\\n")
    script.write("> " + outputDirectory + "/" + sample + ".bed")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
