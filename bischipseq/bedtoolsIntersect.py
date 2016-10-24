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
parser = argparse.ArgumentParser(description='Bedtools intersect narrowPeaks with methylation calls.')
parser.add_argument("-c", "--scriptsDirectory", help="Scripts directory.", default="bedtoolsIntersect")
parser.add_argument("-m", "--macs2Directory", help="MACS2 directory with narrowPeaks files.", default="../macs2/without_duplicates/q_value_20_percent/")
parser.add_argument("-b", "--bismarkDirectory", help="Bismark directory with cov files.", default="../bismark/")
parser.add_argument("-o", "--outputDirectory", help="Output directory with sorted bedgraph files.", default="../bedtoolsIntersect/narrowPeaks_and_cov/")
parser.add_argument("-s", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
macs2Directory = os.path.abspath(args.macs2Directory)
bismarkDirectory = os.path.abspath(args.bismarkDirectory)
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
    scriptName = 'bedtoolsIntersect_' + sample + '.sh'
    script = open(scriptName, 'w')
    util.writeHeader(script, config, "bedtoolsIntersect")
    script.write("bedtools intersect " + "\\\n")
    script.write("-wa -wb " + "\\\n")
    script.write("-a " + macs2Directory +  "/" + sample + "/" + sample + "_peaks.narrowPeak " + "\\\n")
    script.write("-b " + bismarkDirectory + "/" + sample + "/methylation_extractor/" + sample + "_deduplicated_sorted_by_read_name_zero_start.bismark.cov " + "\\\n")
    script.write("> " + outputDirectory + "/" + sample + ".bed")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
