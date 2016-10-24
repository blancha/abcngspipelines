#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates star scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=star", default="star")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files. DEFAULT=../data/FASTQ_files/untrimmed", default="../data/FASTQ_files/untrimmed")
parser.add_argument("-o", "--outputDirectory", help="Output directory with star results. DEFAULT=../results/star", default="../results/star")
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
ppn = config.get("star", "ppn")

genome = config.get("project", "genome")
readLength = config.get("star", "readLength")
starIndex = config.get(genome, "starIndex_readlength" + readLength)
sjdbGTFfile = config.getboolean("star", "sjdbGTFfile")
if sjdbGTFfile:
   sjdbGTFfile = config.get(genome, "gtfFile")
else:
   sjdbGTFfile = "None"
if config.has_option("star", "chimSegmentMin"):
    chimSegment = config.get("star", "chimSegmentMin")
else: 
    chimSegmentMin = None
quantMode = config.get("star", "quantMode")
readFilesCommand = config.get("star", "readFilesCommand")
runThreadN = config.get("star", "runThreadN")
outSAMtype = config.get("star", "outSAMtype")
outWigType = config.get("star", "outWigType")
outWigStrand = config.get("star", "outWigStrand")
outFilterIntronMotifs = config.get("star", "outFilterIntronMotifs")
alignEndsType = config.get("star", "alignEndsType")

# Read samples file.
samplesFile = util.readSamplesFile()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Cycle through all the samples and write the star scripts.
for index, row in samplesFile.iterrows():
    sample = row["sample"]
    file_R1 = row["file_r1"]
    file_R2 = row["file_r2"]
    # Create output directories
    if not os.path.exists(outputDirectory + "/" +  sample):
        os.mkdir(outputDirectory + "/" + sample)    
    # Create script file.
    scriptName = "star_" + sample + ".sh"
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "star")
    script.write("STAR " + "\\\n")
    script.write("--runMode alignReads" + " \\\n")
    script.write("--runThreadN " + runThreadN + " \\\n")
    script.write("--genomeDir " + starIndex + " \\\n")
    script.write("--sjdbOverhang " + str(int(readLength) - 1) + " \\\n")
    if not sjdbGTFfile == "None":
        script.write("--sjdbGTFfile " + sjdbGTFfile + " \\\n")
    script.write("--readFilesIn " + "\\\n")
    script.write(os.path.relpath(os.path.join(inputDirectory, file_R1)) + " \\\n")
    script.write(os.path.relpath(os.path.join(inputDirectory, file_R2)) + " \\\n")
    if readFilesCommand == "zcat" or readFilesCommand == "gunzip -c" or readFilesCommand == "bunzip -c":
               script.write("--readFilesCommand " + readFilesCommand + " \\\n")
    script.write("--outFileNamePrefix " + os.path.relpath(os.path.join(outputDirectory, sample, sample)) + " \\\n")
    if not outSAMtype == "None":
        script.write("--outSAMtype " + outSAMtype + " \\\n")
        if outSAMtype == "BAM SortedByCoordinate":
            limit = str(2700000000 * int(ppn))
            script.write("--limitBAMsortRAM " + limit + " \\\n")
    if not outWigType == "None":
        script.write("--outWigType " + outWigType + " \\\n")
        script.write("--outWigStrand " + outWigStrand + " \\\n")
    if chimSegmentMin is not None:
        script.write("--chimSegmentMin " + chimSegmentMin + "\\\n")
    if not quantMode == "None":
        script.write("--quantMode " + quantMode + " \\\n")
    if not outFilterIntronMotifs == "None":
        script.write("--outFilterIntronMotifs " + outFilterIntronMotifs + " \\\n")
    if not alignEndsType == "None":
        script.write("--alignEndsType " + alignEndsType + " \\\n")

    script.write("&> " + scriptName + ".log")    

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
