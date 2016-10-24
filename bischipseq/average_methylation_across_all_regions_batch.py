#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 04/06/2014

import argparse
import glob
import os
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description='Generates scripts to calculate average methylation levels across regions.')
parser.add_argument("-c", "--scriptsDirectory", help="Scripts directory.", default="average_methylation_across_all_regions")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BED files, containing methylation percentages.", default="../bedtoolsGroupBy/narrowPeaks_and_cov/")
parser.add_argument("-o", "--outputDirectory", help="Output directory with sorted bedgraph files.", default="../average_methylation_across_all_regions")
parser.add_argument("-s", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Read configuration files.
config = util.readConfigurationFiles()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.mkdir(outputDirectory)

for file in os.listdir(inputDirectory):
    file = os.path.splitext(file)[0]
    # Create script file.
    scriptName = 'average_methylation_across_all_regions_' + file + '.sh'
    script = open(scriptName, 'w')
    util.writeHeader(script, config, "average_methylation_across_all_regions")
    script.write("average_methylation_across_all_regions.py " + "\\\n")
    script.write("--input_file " + inputDirectory +  "/" + file + ".bed " + "\\\n")
    script.write("--output_file " + outputDirectory + "/" + file + ".txt")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
