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
parser = argparse.ArgumentParser(description="Generates Picard tools CollectAlignmentSummaryMetrics scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=collectalignmentsummarymetrics", default="collectalignmentsummarymetrics")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/bamtools_filter", default="../results/bamtools_filter")
parser.add_argument("-o", "--outputDirectory", help="Output directory with htseqcount results. DEFAULT=../results/collectalignmentsummarymetrics", default="../results/collectalignmentsummarymetrics")
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

genome = config.get("project", "genome")
genomeFile = config.get(genome, "genomeFile")
picard_folder = config.get("picard", "folder")

samples = util.getMergedsamples()
#samples = glob.glob("DNA*") 

# Create script and output directories, if they do not exist yet.
util.makeDirectory(outputDirectory)
util.makeDirectory(scriptsDirectory)

# CD to scripts directory
os.chdir(scriptsDirectory)

# Write scripts
for sample in samples:
    
    scriptName =  "collectalignmentsummarymetrics_" + sample + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "collectalignmentsummarymetrics")
    # Reorder 
    script.write("java -jar " + os.path.join(picard_folder, "CollectAlignmentSummaryMetrics.jar") + " \\\n")
    script.write("METRIC_ACCUMULATION_LEVEL=ALL_READS" + " \\\n")
    script.write("INPUT=" + os.path.join(inputDirectory, sample + ".filtered.bam") + " \\\n")
    script.write("OUTPUT=" + os.path.join(outputDirectory, sample + "_picard_alignment_metrics.txt") + " \\\n")    
    script.write("REFERENCE_SEQUENCE=" + genomeFile + " \\\n")
    script.write("VALIDATION_STRINGENCY=LENIENT " + "\\\n")
    script.write("&> " + scriptName + ".log")
    
    script.close()
