#!/usr/bin/env python

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 24/03/2014

import argparse
import ConfigParser
import glob
import os
import subprocess
import util

parser = argparse.ArgumentParser(description='Converts BAM files to bedGraph files.')
parser.add_argument("-s", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="yes")
args = parser.parse_args()

# Read configuration files
config = util.readConfigurationFiles()

# Create script directory, if it does not exist yet, and cd to it.
if not os.path.exists("coverage_bedgraph"):
    os.mkdir("coverage_bedgraph")
    
os.chdir("coverage_bedgraph")

# Read samples files
samplesFile = open("../samples.txt")
samplesLines = samplesFile.readlines()[1:]

samples = []

for line in samplesLines:
    samples.append(line.split()[3])
 
# Keep only unique samples.   
samples = list(set(samples))

for sample in samples:
    # Open coverage file.
    coverageFile = glob.glob("../../bismark/" + sample + "/methylation_extractor/" + sample + "_bismark_pe.deduplicated.bismark.cov")[0]
    # Create script file.
    scriptName = 'coverage_bedgraph_' + sample + '.sh'
    script = open(scriptName, 'w')
    util.write_header(script, config, "coverage_bedgraph")
    script.write("coverage_bedgraph.py --file=" + coverageFile)
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
