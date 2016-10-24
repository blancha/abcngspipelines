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
parser = argparse.ArgumentParser(description="Converts BAM files to bedgraph files.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=bamtobedgraph_stranded", default="bamtobedgraph_stranded")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/bamfiles_stranded", default="../results/bamfiles_stranded")
parser.add_argument("-o", "--outputDirectory", help="Output directory with bedgraph files. DEFAULT=../results/bedgraph_stranded", default="../results/bedgraph_stranded")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

header = config.getboolean("server", "PBS_header")
genome = config.get("project", "genome")
chromInfo = config.get(genome, "chromInfo")
UCSCchromInfo = config.get(genome, "UCSCchromInfo")
institute = config.get(genome, "institute")

samples = util.getsamples(lanes=False)

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Output to institute subdirectory
outputDirectory_institute = os.path.join(outputDirectory, institute)
if not os.path.exists(outputDirectory_institute):
        os.makedirs(outputDirectory_institute)

############################
# bamtobedgraph.sh scripts #
############################
for sample in samples:
    scriptName = "bamtobedgraph_" + sample + "_negative.sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "bamtobedgraph")
    # BAM to bedgraph
    script.write("bedtools genomecov" + " \\\n")
    # Only if type == rnaseq ???
    script.write("-split" + " \\\n")
    script.write("-bg" + " \\\n")
    script.write("-g " + chromInfo + " \\\n")
    script.write("-ibam " + os.path.relpath(os.path.join(inputDirectory, sample + "_negative.bam")) + " \\\n")
    script.write("1> " + os.path.relpath(os.path.join(outputDirectory_institute, sample + "_negative.bedgraph")) + " \\\n")
    script.write("2> " + scriptName + ".log")
    script.close()
    scriptName = "bamtobedgraph_" + sample + "_positive.sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "bamtobedgraph")
    # BAM to bedgraph
    script.write("bedtools genomecov" + " \\\n")
    # Only if type == rnaseq ???
    script.write("-split" + " \\\n")
    script.write("-bg" + " \\\n")
    script.write("-g " + chromInfo + " \\\n")
    script.write("-ibam " + os.path.relpath(os.path.join(inputDirectory, sample + "_positive.bam")) + " \\\n")
    script.write("1> " + os.path.relpath(os.path.join(outputDirectory_institute, sample + "_positive.bedgraph")) + " \\\n")
    script.write("2> " + scriptName + ".log")
    script.close()
if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
