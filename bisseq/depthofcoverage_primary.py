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
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=depthofcoverage_primary", default="depthofcoverage_primary")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/bamtools_filter", default="../results/bamtools_filter")
parser.add_argument("-o", "--outputDirectory", help="Output directory with filterd BAM files. DEFAULT=../results/depthofcoverage_primary", default="../results/depthofcoverage_primary")
parser.add_argument("-r", "--target_BED", help="Target BED file. DEFAULT=/sb/project/afb-431/BIF/SZY/SZY-BIF-P3/data/140501_RN4_PTSD3_RM_EPI_combined_coverage/140501_RN4_PTSD3_RM_EPI_primary_targets.bed", default="/sb/project/afb-431/BIF/SZY/SZY-BIF-P3/data/140501_RN4_PTSD3_RM_EPI_combined_coverage/140501_RN4_PTSD3_RM_EPI_primary_targets.bed")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the scripts directory, cd to the scripts directory
util.cdMainScriptsDirectory()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
target_BED = args.target_BED

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

genome = config.get("project", "genome")
genomeFile = config.get(genome, "genomeFile")

samples = util.getMergedsamples()

# Create script and output directories, if they do not exist yet.
util.makeDirectory(outputDirectory)
util.makeDirectory(scriptsDirectory)

# CD to scripts directory
os.chdir(scriptsDirectory)

# Write scripts
for sample in samples:
    scriptName =  "depthofcoverage_primary_" + sample + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "depthofcoverage_primary")
    # Reorder 
    
    script.write("java -Xmx4g -Xms4g -jar /software/areas/ircm/tools/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar" + " \\\n")
    script.write("-T DepthOfCoverage" + " \\\n")
    script.write("-R " + genomeFile + " \\\n")
    script.write("-I " + os.path.join(inputDirectory, sample + ".filtered.bam") + " \\\n")
    script.write("-o " + os.path.join(outputDirectory, sample + "_gatk_primary_target_coverage.txt") + " \\\n")
    script.write("-L  " + target_BED + " \\\n")
    script.write("-ct 1 -ct 10 -ct 20" + " \\\n")
    script.write("&> " + scriptName + ".log")

    script.close()
