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
parser = argparse.ArgumentParser(description="Generates Picard target list scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=collectinsertsizemetrics", default="collectinsertsizemetrics")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/bamtools_merge", default="../results/star")
parser.add_argument("-o", "--outputDirectory", help="Output directory with filterd BAM files. DEFAULT=../results/collectinsertsizemetrics", default="../results/collectinsertsizemetrics")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the scripts directory, cd to the scripts directory
util.cdMainScriptsDirectory()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
scriptsDirectory = os.path.abspath(args.scriptsDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

picard_folder = config.get("picard", "folder")
genome = config.get("project", "genome")
genomeFile = config.get(genome, "genomeFile")

# Read samples file.
samplesFile = util.readSamplesFile()
samples = samplesFile["sample"]

# Create script and output directories, if they do not exist yet.
util.makeDirectory(outputDirectory)
util.makeDirectory(scriptsDirectory)

# CD to scripts directory
os.chdir(scriptsDirectory)

# Write scripts
for sample in samples:
    scriptName =  "collectinsertsizemetrics_" + sample + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "collectinsertsizemetrics")
    # Reorder 
    
    script.write("java -Xmx4g -jar " + os.path.join(picard_folder, "CollectInsertSizeMetrics.jar") + " \\\n")
    script.write("VALIDATION_STRINGENCY=LENIENT " + "\\\n")
    script.write("HISTOGRAM_FILE=" + os.path.join(outputDirectory, sample + "_picard_insert_size_plot.pdf") + " \\\n")
    script.write("INPUT=" + os.path.join(inputDirectory, sample + ".filtered.bam") + " \\\n")
    script.write("OUTPUT=" + os.path.join(outputDirectory, sample + "_picard_insert_size_metrics.txt") + " \\\n")
    script.write("&> " + scriptName + ".log")

    script.close()
