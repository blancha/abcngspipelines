#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 21/05/2014

import argparse
import glob
import os
import subprocess
import util

parser = argparse.ArgumentParser(description='Generates script to extract methylation counts from Bismark coverage file, and create methylated and unmethylated bedgraph files.')
parser.add_argument("-s", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

samples = util.getMergedsamples()

# Read configuration files
config = util.readConfigurationFiles()

# Create script directory, if it does not exist yet, and cd to it.
if not os.path.exists("extract_methylation_counts_bismark_cov"):
    os.mkdir("extract_methylation_counts_bismark_cov")
    
os.chdir("extract_methylation_counts_bismark_cov")

# Create bedgraph directory, if it does not exist yet.
outputDirectory = "../../bedgraph/methylation_counts/"
if not os.path.exists(outputDirectory):
    os.mkdir(outputDirectory)

inputDirectory = "../../bismark/"

for sample in samples:
    coverage_file = glob.glob(inputDirectory + sample + "/methylation_extractor/*.cov")[0]
    methylated_file = outputDirectory + sample + "_methylated_counts.bedgraph"
    unmethylated_file = outputDirectory + sample + "_unmethylated_counts.bedgraph"
    # Create script file.
    scriptName = 'extract_methylation_counts_bismark_cov_' + sample + '.sh'
    script = open(scriptName, 'w')
    util.writeHeader(script, config, "extract_methylation_counts_bismark_cov")
    script.write("extract_methylation_counts_bismark_cov.py " + "\\\n")
    script.write("--coverage " + coverage_file + " \\\n")
    script.write("--methylated " + methylated_file + " \\\n")
    script.write("--unmethylated " + unmethylated_file)
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
