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
parser = argparse.ArgumentParser(description="Generates Trimmomatic scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="trimmomatic")
parser.add_argument("-o", "--outputDirectory", help="Output directory with trimmed FASTQ files.", default="../data/FASTQ_files/trimmed")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
outputDirectory = os.path.abspath(args.outputDirectory)

# Read samples file.
samplesDataFrame = util.readSamplesFile()
samples = samplesDataFrame["sample"].tolist()

# Create samples file after trimming
samples_after_trimming = open("samples_after_trimming.txt", "w")
samples_after_trimming.write("\t".join(samplesDataFrame.columns.tolist()))
samples_after_trimming.write("\n")

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
trimAdapters = config.getboolean("trimmomatic", "trimAdapters")
toolsFolder = config.get("server", "toolsFolder")
threads = config.get("trimmomatic", "threads")
minlength = config.get("trimmomatic", "minlength")
programFile = config.get("trimmomatic", "programFile")
adaptersFile = config.get("trimmomatic", "adaptersFile")

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory + "/unpaired"):
    os.makedirs(outputDirectory + "/unpaired")

# Cycle through all the samples and write the fastqc scripts.
for index, row in samplesDataFrame.iterrows():
    sample = row["sample"]
    file_r1 = row["file_r1"]
    file_r2 = row["file_r2"]
    group = row["group"]
    samples_after_trimming.write(
        os.path.relpath(os.path.join(outputDirectory, os.path.basename(file_r1))) + "\t" +
        os.path.relpath(os.path.join(outputDirectory, os.path.basename(file_r2))) + "\t" +
        sample + "\t" + group + "\n")
    files = [file_r1, file_r2]
    if "lane" in samplesDataFrame.columns:
        lane = "_" + row["lane"]
    else:
        lane = ""
    sample +=lane
    # Create script file.
    scriptName = 'trimmomatic_' + sample + '.sh'
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "trimmomatic")
    script.write("java -jar " + programFile + " \\\n")
    script.write("PE" + " \\\n")
    script.write("-threads " + threads + " \\\n")
    script.write(os.path.relpath(file_r1) + " \\\n")
    script.write(os.path.relpath(file_r2) + " \\\n")
    script.write(os.path.relpath(os.path.join(outputDirectory, os.path.splitext(os.path.basename(file_r1))[0] + lane + os.path.splitext(os.path.basename(file_r1))[1])) + " \\\n")
    script.write(os.path.relpath(os.path.join(outputDirectory, "unpaired", os.path.splitext(os.path.basename(file_r1))[0] + lane + os.path.splitext(os.path.basename(file_r1))[1])) + " \\\n")
    script.write(os.path.relpath(os.path.join(outputDirectory, os.path.splitext(os.path.basename(file_r2))[0] + lane + os.path.splitext(os.path.basename(file_r2))[1])) + " \\\n")
    script.write(os.path.relpath(os.path.join(outputDirectory, "unpaired", os.path.splitext(os.path.basename(file_r2))[0] + lane + os.path.splitext(os.path.basename(file_r2))[1])) + " \\\n")
    if trimAdapters:
        script.write("ILLUMINACLIP:" + adaptersFile + ":2:30:10 " + "\\\n")
    script.write("LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:" + minlength + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.close()    

samples_after_trimming.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
