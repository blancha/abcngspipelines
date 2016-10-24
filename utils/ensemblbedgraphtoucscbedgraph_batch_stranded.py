#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: ../2014

import argparse
import glob
import os
import pandas
import subprocess
import util

parser = argparse.ArgumentParser(description='Generate scripts to convert bedgraph files from Ensembl to UCSC format.')
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=ensemblbedgraphtoucscbedgraph", default="ensemblbedgraphtoucscbedgraph")
parser.add_argument("-i", "--inputDirectory", help="Input directory with Ensembl bedgraph files. DEFAULT=../results/bedgraph/Ensembl/", default="../results/bedgraph/Ensembl/")
parser.add_argument("-o", "--outputDirectory", help="Output directory with UCSC bedgraph files. DEFAULT=../results/bedgraph/UCSC/", default="../results/bedgraph/UCSC/")
parser.add_argument("-d", "--dictionary", help="Tab-delimited text file with the Ensembl symbols and the corresponding UCSC symbols. Default=/sb/project/afb-431-ac/genomes/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/ucsctoEnsembl.txt", default="/sb/project/afb-431-ac/genomes/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/ucsc_to_ensembl_chromsomes.dict")
parser.add_argument("-t", "--stranded", help="Stranded bedgraph files. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no") 
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory)
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
if (args.stranded == "yes") | (args.stranded == "y"):
    stranded = True
else:
    stranded = False
dictionary = os.path.abspath(args.dictionary)

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

genome = config.get("project", "genome")

samplesFile = pandas.read_csv("samples.txt", sep="\t")

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Cycle through all the samples and write the scripts.
for index, row in samplesFile.iterrows():
    sample = row["sample"]
    # Create script file.
    scriptName = 'ensemblbedgraphtoucscbedgraph_' + sample + '.sh'
    script = open(scriptName, 'w')
    util.writeHeader(script, config, "ensemblbedgraphtoucscbedgraph")
    script.write("ensemblbedgraphtoucscbedgraph.py " + "\\\n")
    script.write("--ensembl_bedgraph " + os.path.join(inputDirectory,sample + ".bedgraph") + " \\\n")
    script.write("--ucsc_bedgraph " + os.path.join(outputDirectory, sample + ".bedgraph") + " \\\n")
    script.write("--dictionary " + dictionary + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.close()
    if stranded: 
        # Create positive strand script file.
        scriptName = 'ensemblbedgraphtoucscbedgraph_' + sample + '_positive.sh'
        script = open(scriptName, 'w')
        util.writeHeader(script, config, "ensemblbedgraphtoucscbedgraph")
        script.write("ensemblbedgraphtoucscbedgraph.py " + "\\\n")
        script.write("--ensembl_bedgraph " + os.path.join(inputDirectory,sample + "_positive.bedgraph") + " \\\n")
        script.write("--ucsc_bedgraph " + os.path.join(outputDirectory, sample + "_positive.bedgraph") + " \\\n")
        script.write("--dictionary " + dictionary + " \\\n")
        script.write("&> " + scriptName + ".log")
        script.close()
        # Create negative strand script file.
        scriptName = 'ensemblbedgraphtoucscbedgraph_' + sample + '_negative.sh'
        script = open(scriptName, 'w')
        util.writeHeader(script, config, "ensemblbedgraphtoucscbedgraph")
        script.write("ensemblbedgraphtoucscbedgraph.py " + "\\\n")
        script.write("--ensembl_bedgraph " + os.path.join(inputDirectory,sample + "_negative.bedgraph") + " \\\n")
        script.write("--ucsc_bedgraph " + os.path.join(outputDirectory, sample + "_negative.bedgraph") + " \\\n")
        script.write("--dictionary " + dictionary + " \\\n")
        script.write("&> " + scriptName + ".log")
        script.close()
if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
