#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 12/04/2014

import argparse
import glob
import os
import os.path
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates bamtofastq scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="bamtofastq")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../data/BAM_files", default="../data/BAM_files")
parser.add_argument("-o", "--outputDirectory", help="Output directory with FASTQ files DEFAULT=../data/FASTQ_files/untrimmed.", default="../data/FASTQ_files/untrimmed")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(inputDirectory)

# Create script and output directories, if they do not exist yet.
util.makeDirectory(outputDirectory)
util.makeDirectory(scriptsDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
picard_folder = config.get("picard", "folder")

# Change to scripts directory
os.chdir(scriptsDirectory)

# Store all the BAM filenames.
files = glob.glob(inputDirectory + "/*bam")

# Cycle through all the files, creating the scripts.
for file in files:
    sample = os.path.basename(file)[:-4]
    scriptName = 'bamtofastq_' + sample + '.sh'
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "bamtofastq")

    script.write("java -Xmx4g -Xms4g -jar " + os.path.join(picard_folder, "picard.jar") + " SamToFastq" + " \\\n")
    script.write("VALIDATION_STRINGENCY=SILENT " + "\\\n")
    script.write("INPUT=" + file + " \\\n")
    script.write("FASTQ=" + os.path.join(outputDirectory,  sample + "_R1.fastq") + " \\\n")
    script.write("SECOND_END_FASTQ=" + os.path.join(outputDirectory, sample + "_R2.fastq")  + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.write("\n\n")

    script.write("gzip " + os.path.join(outputDirectory, sample + "_R1.fastq") + " \\\n")
    script.write("&>> " + scriptName + ".log")
    script.write("\n\n")
    script.write("gzip " + os.path.join(outputDirectory, sample + "_R2.fastq") + " \\\n")
    script.write("&>> " + scriptName + ".log")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
