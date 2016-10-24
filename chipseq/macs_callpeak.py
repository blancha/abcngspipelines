#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
from collections import OrderedDict
import glob
import os
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates macs_callpeak scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="macs_callpeak")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files.", default="../results/bowtie/")
parser.add_argument("-o", "--outputDirectory", help="Output directory with macs callpeak results.", default="../results/macs_callpeak/")
parser.add_argument("-u", "--subDirectories", help="BAM files located in sample subdirectories in inputDirectory. DEFAULT=yes", choices=["yes", "no", "y", "n"], default="yes")
parser.add_argument("-p", "--pipeline", help="Write script as part of pipeline, with files in known location")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
pipeline = args.pipeline
subdirectories = args.subDirectories

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")

genome = config.get("project", "genome")
# Get all the comparisons.
if config.has_option("macs_callpeak", "comparisons"):
    comparisons = config.get("macs_callpeak", "comparisons").split(",")
else:
    print("Error!\nThe configuration has no macs section, with a comparisons subsection.")
    exit(1)
if config.has_option("macs_callpeak", "broad"):
    broad = config.getboolean("macs_callpeak", "broad")
else:
    print("Error!\nThe configuration has no macs_callpeak section, with a broad subsection.")
    exit(1)
if config.has_option("macs_callpeak", "keep-dup"):
    keepdup = config.get("macs_callpeak", "keep-dup")
else:
    print("Error!\nThe configuration has no macs_callpeak section, with a keep-dup subsection.")
    exit(1)

qvalue = config.get('macs_callpeak', 'qvalue')

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)


# Get the genome size
if (genome == "mm9") | (genome=="mm10") | (genome == "GRCm38") | (genome == "NCBIM37"):
    genomeSize = "mm"
elif (genome == "hg19") | (genome == "GRCh37"):
    genomeSize = "hs"
elif (genome == "dm3") | (genome == "dm6") | (genome == "BDGP5") | (genome == "BDGP6"):
    genomeSize = "dm"

# Cycle through all the comparisons and write the macs_callpeak scripts.
for comparison in comparisons:
    samples = comparison.split("vs")
    treatment = samples[0].strip()
    if len(samples) == 2:
        control = samples[1].strip()
        comparisonName = treatment + "_vs_" + control
    else:
        comparisonName = treatment
    # Create output directory for each comparison
    if not os.path.exists(outputDirectory + "/" + comparisonName):
        os.makedirs(outputDirectory + "/" + comparisonName)
    if pipeline:
        treatmentFile = os.path.join(inputDirectory, treatment, treatment + ".bam")
        if(len(samples) == 2):
            controlFile = os.path.join(inputDirectory, control, control + ".bam")
    else:
        if (subdirectories == "yes") | (subdirectories=="y"):
            treatmentFile = os.path.join(inputDirectory, treatment, treatment + ".bam")
        else:
            treatmentFile = os.path.join(inputDirectory, treatment + ".bam")
        if len(samples) == 2:
            if (subdirectories == "yes") | (subdirectories=="y"):
                controlFile = os.path.join(inputDirectory, control, control + ".bam")
            else:
                controlFile = os.path.join(inputDirectory, control + ".bam")
    # Create script file.
    scriptName = 'macs_callpeak_' + comparisonName + '.sh'
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "macs_callpeak")
    #Deactivate Python 3 virtual environment, and activate Python 2 virtual environment to be able to use TopHat.
    script.write("source /sb/software/areas/ircm/tools/python_environments/python2.7/bin/activate\n\n")
    script.write("macs2 callpeak " + "\\\n")
    if broad:   
        script.write("--broad " + "\\\n")
    script.write("--bdg " + "\\\n")
    script.write("--gsize " + genomeSize + " \\\n")
    script.write("--keep-dup " + keepdup + " \\\n")
    script.write("--qvalue " + qvalue + " \\\n") 
    script.write("--name " + outputDirectory + "/" + comparisonName + "/" + comparisonName + " " + "\\\n")
    script.write("--treatment " + treatmentFile + " \\\n")
    if len(samples) == 2:
        script.write("--control " + controlFile + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)

