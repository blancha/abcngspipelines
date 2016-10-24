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
parser = argparse.ArgumentParser(description="Generates Cutadapt scripts for small RNA-Seq data.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=cutadapt", default="cutadapt")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files. DEFAULT=../data/FASTQ_files/untrimmed", default="../data/FASTQ_files/untrimmed")
parser.add_argument("-o", "--outputDirectory", help="Output directory with trimmed FASTQ files. DEFAULT=../data/FASTQ_files/trimmed", default="../data/FASTQ_files/trimmed")
parser.add_argument("-z", "--gzip", help="gzip output folders. DEFAULT=yes", choices=["yes", "no", "y", "n"], default="yes")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
adapter = config.get("cutadapt", "adapter")
minlength = config.get("cutadapt", "minimum-length")
maxlength = config.get("cutadapt", "maximum-length")
qualitycutoff = config.get("cutadapt", "quality-cutoff")
trimn = config.getboolean("cutadapt", "trim-n")
cut = config.get("cutadapt", "cut")
tooshortoutput = config.getboolean("cutadapt", "tooshortoutput")
toolongoutput = config.getboolean("cutadapt", "toolongoutput")

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)
    os.makedirs(os.path.join(outputDirectory, "too_short_after_trimming"))
    os.makedirs(os.path.join(outputDirectory, "too_long_after_trimming"))

# Store the list of files with the extensions fastq or fastq.gz 
files = glob.glob(inputDirectory + "/*.fastq") + glob.glob(inputDirectory + "/*.fastq.gz")
files.sort()

# Write the script(s)
# Cycle through all the R1 files.
for file in files: 
    fileR1=os.path.basename(file)
    # Create script file.
    scriptName = 'cutadapt_' + fileR1 + '.sh'
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "cutadapt")
    script.write("cutadapt" + " \\\n")
    if not adapter == "None":
        script.write("--adapter " + adapter + " \\\n")
    if not minlength == "None":
        script.write("--minimum-length " + minlength + " \\\n")
    if not maxlength == "None":
        script.write("--maximum-length " + maxlength + " \\\n")
    if not qualitycutoff == "None":
        script.write("--quality-cutoff " + qualitycutoff + " \\\n")
    if trimn:
        script.write("--trim-n" + " \\\n")
    if not cut == "None":
            script.write("--cut " + cut + " \\\n")
    if (args.gzip.lower() == "no") | (args.gzip.lower() == "n"):
        script.write("--output " + os.path.relpath(os.path.join(outputDirectory, fileR1[:-3])) + " \\\n")
    else: 
        script.write("--output " + os.path.relpath(os.path.join(outputDirectory, fileR1)) + " \\\n")
    if tooshortoutput:
        script.write("--too-short-output " + os.path.relpath(os.path.join(outputDirectory, "too_short_after_trimming", fileR1)) + " \\\n")
    if toolongoutput:
        script.write("--too-long-output " + os.path.relpath(os.path.join(outputDirectory, "too_long_after_trimming", fileR1)) + " \\\n")
    script.write(os.path.relpath(os.path.join(inputDirectory, fileR1)) + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.close()    

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
