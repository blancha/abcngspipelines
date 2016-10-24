#!/usr/bin/env python

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 10/03/2014

import argparse
import ConfigParser
import glob
import os
import subprocess
import sys
sys.path.append('/sb/software/areas/ircm/tools/ngs_pipelines/utils/')
import util

parser = argparse.ArgumentParser(description='Generate methylation_extractor scripts.')
parser.add_argument("-s", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="yes")
args = parser.parse_args()

# Read configuration file
config = util.readConfigurationFiles()

# Create script directory, if it does not exist yet, and cd to it.
if not os.path.exists("methylation_extractor"):
    os.mkdir("methylation_extractor")
    
os.chdir("methylation_extractor")

# Cycle through all the samples and write the methylation_extractor scripts.
directories = os.listdir("../../bismark")

samplesFile = open("../samples.txt")
samplesLine = samplesFile.readlines()[1::2]

samples = []

for line in samplesLine:
    samples.append(line.split()[3].split("_")[-2])
    
uniquesamples = list(set(samples))


for sample in uniquesamples:
    # Create output directory
    outputDirectory = "../../bismark/" + sample + "/methylation_extractor"
    if not os.path.exists(outputDirectory):
        os.mkdir(outputDirectory)
    bamFile = str(glob.glob("../../bismark/" + sample + "/*.bam"))[2:-2]
    # Create script file.
    scriptName = 'methylation_extractor_' + sample + '.sh'
    script = open(scriptName, 'w')
    util.write_header(script, config, "methylation_extractor")
    script.write("bismark_methylation_extractor --buffer_size 10G --paired --no_overlap --comprehensive --report " + "\\\n")
    script.write("--bedGraph --CX --cutoff " + config.get('methylation_extractor', 'cutoff') + " " + "\\\n")
    script.write("--output " + outputDirectory + " " + "\\\n")
    script.write(bamFile)
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)

