#!/usr/bin/env python3


# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates scripts to download SRA files.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=.", default=".")
parser.add_argument("-i", "--samplesFile", help="Input file with names of SRA runs. DEFAULT=.", default="./SraRunTable.txt")
parser.add_argument("-o", "--outputDirectory", help="Output directory with SRA files. DEFAULT=.", default=".")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
#util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
samplesFile = os.path.abspath(args.samplesFile)
outputDirectory = os.path.abspath(args.outputDirectory)

# Check if the samplesFile exists, and is a file.
if not(os.path.exists(samplesFile) and os.path.isfile(samplesFile)):
    exit(samplesFile +  " does not exist or is not a file.")

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

# Read input file.
samplesFile = pandas.read_csv(samplesFile, sep="\t")

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Cycle through all the samples and write the star scripts.
for index, row in samplesFile.iterrows():
    run = row["Run_s"]
    # Create script file.
    scriptName = "getsra_" + run + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "getsra")
    script.write("wget" + " \\\n")
    script.write("ftp://ftp.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/" + os.path.join(run[0:6], run, run + ".sra") + " \\\n")

    script.write("&> " + scriptName + ".log")    

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
