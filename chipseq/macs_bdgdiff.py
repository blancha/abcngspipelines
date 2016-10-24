#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 11/03/2014

import argparse
import glob
import os
import subprocess
import util

parser = argparse.ArgumentParser(description='Generate macs_callpeak scripts.')
parser.add_argument("-s", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="yes")
args = parser.parse_args()

# Read configuration file
config = util.readConfigurationFiles()

header = config.getboolean("server", "header")
footer = config.getboolean("server", "footer")

genome = config.get("project", "genome")
keepdup = config.get('macs_callpeak', 'keep-dup')
qvalue = config.get('macs_callpeak', 'qvalue')
format = config.get('macs_callpeak', 'format')

# Create script directory, if it does not exist yet, and cd to it.
if not os.path.exists("macs_callpeak"):
    os.mkdir("macs_callpeak")
    
os.chdir("macs_callpeak")

# Get the genome size
if (genome == "mm9") | (genome=="mm10") | (genome == "GRCM38") | (genome == "NCBIM37"):
    genomeSize = "mm"
if (genome == "hg19") | (genome == "GRCh37"):
    genomeSize = "hs"

# Cycle through all the comparisons and write the macs_callpeak scripts.
comparisons = config.get("macs_callpeak", "comparisons").split(",")

for comparison in comparisons:
    samples = comparison.split("vs")
    treatment = samples[0].strip()
    if len(samples) == 2:
        control = samples[1].strip()
        comparisonName = treatment + "_vs_" + control
    else:
        comparisonName = treatment
    # Create output directory
    outputDirectory = "../../macs_callpeak/" + comparisonName
    if not os.path.exists(outputDirectory):
        os.makedirs(outputDirectory)
    treatmentFile = glob.glob("../../bowtie2/" + treatment + "/*" + treatment + ".bam")[0]
    if len(samples) == 2:
        controlFile = glob.glob("../../bowtie2/" + control + "/*" + control + ".bam")[0]
    # Create script file.
    scriptName = 'macs_callpeak_' + comparisonName + '.sh'
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "macs_callpeak")
    script.write("macs2 callpeak " + "\\\n")
    script.write("--gsize " + genomeSize + " \\\n")
    script.write("--keep-dup " + keepdup + " \\\n")
    script.write("--qvalue " + qvalue + " \\\n") 
    script.write("--format " + format + " \\\n")
    script.write("--name ../../macs_callpeak/" + comparisonName + "/" + comparisonName + " " + "\\\n")
    script.write("--treatment " + treatmentFile)
    if len(samples) == 2:
        script.write(" " + "\\\n")
        script.write("--control " + controlFile)
    if footer:
        util.writeFooter(script)
    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)

