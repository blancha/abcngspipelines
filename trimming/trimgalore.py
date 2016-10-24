#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates trim_galore scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="trimgalore")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files.", default="../data/FASTQ_files/untrimmed")
parser.add_argument("-o", "--outputDirectory", help="Output directory with trimmed FASTQ files.", default="../data/FASTQ_files/trimmed")
parser.add_argument("-f", "--fastqcoutputDirectory", help="Output directory with FASTQC report.", default="../results/fastqc/after_trimming")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
fastqcoutputDirectory = os.path.abspath(args.fastqcoutputDirectory)

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
trimAdapters = config.getboolean("trimmomatic", "trimAdapters")
toolsFolder = config.get("server", "toolsFolder")
threads = config.get("trimmomatic", "threads")
minlength=config.get("trimmomatic", "minlength")

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directories, if they do not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)
if not os.path.exists(fastqcoutputDirectory):
    os.makedirs(fastqcoutputDirectory)

# Store the list of files with the extensions fastq or fastq.gz 
files = glob.glob(inputDirectory + "/*.fastq") + glob.glob(inputDirectory + "/*.fastq.gz")
files.sort()

# Write the script(s)
# Cycle through all the files, 2 by 2.
for i in range(0, len(files), 2): 
    fileR1=os.path.basename(files[i])
    fileR2=os.path.basename(files[i+1])
    # Create script file.
    scriptName = 'trimgalore_' + fileR1.replace("_R1", "") + '.sh'
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "trimmomatic")
    script.write("trim_galore" + " \\\n")
    script.write("--fastqc" + " \\\n")
    script.write('--fastqc_args "--outdir ' + fastqcoutputDirectory + '"' + " \\\n")
    script.write("--paired" + " \\\n")
    script.write("--length 50" + " \\\n")
    script.write("--output_dir " + outputDirectory + " \\\n")
    script.write(os.path.join(inputDirectory, fileR1) + " \\\n")
    script.write(os.path.join(inputDirectory, fileR2) + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.close()    

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
