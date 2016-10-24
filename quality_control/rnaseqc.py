#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates rnaseqc scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="rnaseqc")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files.", default="../results/tophat")
parser.add_argument("-o", "--outputDirectory", help="Output directory with rnaseqc results.", default="../results/rnaseqc")
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

processors = config.get("rnaseqc", "processors")

trim = config.getboolean("project", "trim")
stranded = config.getboolean("project", "stranded")
genome = config.get("project", "genome")
genomeFile = config.get(genome, "genomeFile")
gtfFile = config.get(genome, "gtfFile")
rnaseqc_gtfFile = config.get(genome, "rnaseqc_gtfFile")
rrna = config.get(genome, "rrna")

samplesFile = pandas.read_csv("samples.txt", sep="\t")

# Get samples
samples = samplesFile["sample"]
conditions = samplesFile["condition"]

# Change to scripts directory
os.chdir(scriptsDirectory)

###########################
# reorder_index_sample.sh #
###########################

# Create scripts subdirectory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory + "/reorder_index"):
    os.mkdir(scriptsDirectory + "/reorder_index")
os.chdir(scriptsDirectory + "/reorder_index")

for sample in samples:
    scriptName =  "reorder_index_" + sample + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "reorder_index")
    # Reorder 
    script.write("java -jar -Xmx" + str(int(processors) * 2700) + "m " + config.get('picard', 'folder') + "ReorderSam.jar \\\n")
    script.write("I=" + inputDirectory + "/" + sample + "/accepted_hits.bam \\\n")
    script.write("OUTPUT=" + os.path.join(inputDirectory, sample, "accepted_hits_reordered.bam") + " \\\n")
    script.write("REFERENCE=" + genomeFile[:-3] + ".reordered.karyotypic.fa" + " \\\n")
    script.write("CREATE_INDEX=true" + " \\\n")
    script.write("&> " + scriptName + ".log")

os.chdir("..")

##############
# rnaseqc.sh #
##############
scriptName = "rnaseqc.sh"
script = open(scriptName, "w")
if header:
   util.writeHeader(script, config, "rnaseqc")
script.write("java -jar -Xmx" + str(int(processors) * 2700) + "m " + config.get('rnaseqc', 'file') + " \\\n")
script.write("-o " + outputDirectory + " \\\n")
script.write("-BWArRNA " + rrna + " \\\n")
script.write("-r " + genomeFile[:-3] + ".reordered.karyotypic.fa" + " \\\n")
script.write("-t " + rnaseqc_gtfFile + " \\\n")
script.write("-s sampleFile.txt " + "\\\n")
script.write("&> " + scriptName + ".log")
script.write("\n\n")
script.write("rm " + os.path.join(inputDirectory, "*", "*reordered.bam") + " \\\n")
script.write("&>> " + scriptName + ".log")
script.write("\n")
script.write("rm " + os.path.join(inputDirectory, "*", "*reordered.bai") + " \\\n")
script.write("&>> " + scriptName + ".log")
script.close()

###################
# sampleFile.txt. #
###################
sampleFile = open("sampleFile.txt", "w")
sampleFile.write("sample ID\tBam File\tNotes\n")
for index, sample in enumerate(samples):
    sampleFile.write(sample + "\t" + os.path.join(inputDirectory, sample, "accepted_hits_reordered.bam") + "\t" + conditions[index] + "\n")
sampleFile.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
