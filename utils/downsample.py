#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 21/05/2014

import argparse
import glob
import os
import subprocess
import util

# Get the command line arguments
parser = argparse.ArgumentParser(description='Down sample bigWig files with windows of 10 bases.')
parser.add_argument("-s", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# Read configuration files
config = util.readConfigurationFiles()

scriptsDirectory = "downsample"
inputDirectory = "../../bigwig/with_normalisation/all_positions/"
outputDirectory = "../../bigwig/with_normalisation/window_10_bases/"

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)  
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.mkdir(outputDirectory)

for file in os.listdir(inputDirectory):
    # Create script file.
    scriptName = 'downsample_' + file + '.sh'
    script = open(scriptName, 'w')
    util.writeHeader(script, config, scriptsDirectory)
    script.write("toolRunner.sh wigmath.Downsample " + "\\\n")
    script.write("--input " + inputDirectory + file + " \\\n")
    script.write("--output " + outputDirectory + file + " \\\n")
    script.write("--window 10" + "\n")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
