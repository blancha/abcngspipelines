#!/usr/bin/env python3

# Version 1.1
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
parser = argparse.ArgumentParser(description="Generates GATK BaseRecalibrator scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=baserecalibrator", default="baserecalibrator")
parser.add_argument("-i", "--inputDirectory", help="Input directory with BAM files. DEFAULT=../results/bwa", default="../results/bwa")
parser.add_argument("-o", "--outputDirectory", help="Output directory with realigned BAM files. DEFAULT=../results/bwa", default="../results/bwa")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
scriptsDirectory = os.path.abspath(args.scriptsDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
toolsFolder = config.get("server", "toolsFolder")

genome = config.get("project", "genome")
genomeFolder = config.get(genome, "genomeFolder")
genomeFile = config.get(genome, "genomeFile")
xmx = config.get("baserecalibrator", "xmx")

# Get samples
samples = util.getsamples(lanes=True)

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Write the scripts
for sample in samples:
    # Write the script
    scriptName =  "baserecalibrator_" + sample + ".sh"
    script = open(scriptName, "w")
    if header:
        util.writeHeader(script, config, "indelrealigner")
    script.write("java -Xmx" + xmx + " \\\n")
    script.write("-jar " + os.path.join(toolsFolder, "GenomeAnalysisTK.jar") + " \\\n")
    script.write("--analysis_type BaseRecalibrator" + " \\\n")
    script.write("--reference_sequence " + genomeFile + " \\\n")
    script.write("--input_file " + os.path.join(inputDirectory, sample, sample + "_realigned_reads.bam") + " \\\n")
    script.write("--knownSites " + os.path.join(genomeFolder, "1000G_phase1.indels.b37.vcf") + " \\\n")    
    script.write("--knownSites " + os.path.join(genomeFolder, "Mills_and_1000G_gold_standard.indels.b37.vcf") + " \\\n")   
    script.write("--out recal_data.table"  + " \\\n") 
    script.write("&> " + scriptName + ".log")
    script.write("\n\n")

    script.write("java -Xmx" + xmx + " \\\n")
    script.write("-jar " + os.path.join(toolsFolder, "GenomeAnalysisTK.jar") + " \\\n")
    script.write("--analysis_type BaseRecalibrator" + " \\\n")
    script.write("--reference_sequence " + genomeFile + " \\\n")
    script.write("--input_file " + os.path.join(inputDirectory, sample, sample + "_realigned_reads.bam") + " \\\n")
    script.write("--knownSites " + os.path.join(genomeFolder, "1000G_phase1.indels.b37.vcf") + " \\\n")    
    script.write("--knownSites " + os.path.join(genomeFolder, "Mills_and_1000G_gold_standard.indels.b37.vcf") + " \\\n")
    script.write("--BQSR recal_data.table" + " \\\n")   
    script.write("--out post_recal_data.table"  + " \\\n") 
    script.write("&>> " + scriptName + ".log")
    script.write("\n\n")

    script.write("java -Xmx" + xmx + " \\\n")
    script.write("-jar " + os.path.join(toolsFolder, "GenomeAnalysisTK.jar") + " \\\n")
    script.write("--analysis_type AnalyzeCovariates" + " \\\n")
    script.write("--reference_sequence " + genomeFile + " \\\n")
    script.write("--beforeReportFile recal_data.table" + " \\\n")
    script.write("--afterReportFile post_recal_data.table" + " \\\n")
    script.write("--plotsReportFile recalibration_plots.pdf" + " \\\n")
    script.write("&>> " + scriptName + ".log")
    script.write("\n\n")

    script.write("java -Xmx" + xmx + " \\\n")
    script.write("-jar " + os.path.join(toolsFolder, "GenomeAnalysisTK.jar") + " \\\n")
    script.write("--analysis_type PrintReads" + " \\\n")
    script.write("--reference_sequence " + genomeFile + " \\\n")
    script.write("--input_file " + os.path.join(inputDirectory, sample, sample + "_realigned_reads.bam") + " \\\n")
    script.write("--BQSR recal_data.table" + " \\\n")
    script.write("--out " + os.path.join(outputDirectory, sample, sample + "_recal_reads.bam") + " \\\n")
    script.write("&>> " + scriptName + ".log")

    script.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
