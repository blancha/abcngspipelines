#!/usr/bin/env python3

import argparse
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates mirdeep scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory. DEFAULT=mirdeep", default="mirdeep")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files. DEFAULT=../data/FASTQ_files/trimmed", default="../data/FASTQ_files/trimmed")
parser.add_argument("-o", "--outputDirectory", help="Output directory with mirdeep scripts. DEFAULT=../results/mirdeep", default="../results/mirdeep")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()
genome = config.get("project", "genome")
genomeFolder = config.get(genome, "genomeFolder")
genomeFile = config.get(genome, "genomeFile")
bowtieIndex = config.get(genome, "bowtie1Index")
header = config.getboolean("server", "PBS_header")

if not os.path.exists(outputDirectory):
    os.mkdir(outputDirectory)

os.chdir(outputDirectory)

samples =  pandas.read_csv("../../scripts/samples.txt", sep="\t")

for index, row in samples.iterrows():
    # Create directories
    if not os.path.exists(row["sample"]):
        os.mkdir(row["sample"])
    os.chdir(row["sample"])
    # Symbolic links to FASTQ files
    subprocess.call("ln -s " + os.path.join(inputDirectory, row["file_r1"][:-3]), shell=True)
    # Mapper script
    mapper_script = open("mapper.sh", "w")
    if header:
        util.writeHeader(mapper_script, config, "mirdeep")
    mapper_script.write("mapper.pl \\\n")
    mapper_script.write("*.fastq \\\n")
    mapper_script.write("-h -n -o 4 -e -m -v \\\n")
    mapper_script.write("-p " + bowtieIndex + " \\\n")
    mapper_script.write("-s " + row["sample"] + "_collapsed.fa \\\n")
    mapper_script.write("-t " + row["sample"] + "_collapsed_vs_genome_" + genome + ".arf" + " \\\n")
    mapper_script.write("1> mapper.sh_output \\\n")
    mapper_script.write("2> mapper.sh_error")
    mapper_script.write("\n\n")
    mapper_script.close()
    if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
        subprocess.call("submitJobs.py", shell=True)
    # Mirdeep script
    mirdeep_script = open("mirdeep.sh", "w")
    if header:
        util.writeHeader(mirdeep_script, config, "mirdeep")
    mirdeep_script.write("miRDeep2.pl \\\n")
    mirdeep_script.write(row["sample"] + "_collapsed.fa \\\n")
    mirdeep_script.write(genomeFile + " \\\n")
    mirdeep_script.write(row["sample"] + "_collapsed_vs_genome_" + genome + ".arf" + " \\\n")
    mirdeep_script.write(os.path.join(genomeFolder, "mirbase21", "mature.hsa.dna.fa") + " \\\n")
    mirdeep_script.write(os.path.join(genomeFolder, "mirbase21", "mature.ptr.dna.fa") + " \\\n")
    mirdeep_script.write(os.path.join(genomeFolder, "mirbase21", "hairpin.hsa.dna.fa") + " \\\n")
    mirdeep_script.write("-t Human \\\n")
    mirdeep_script.write("-P \\\n")
    mirdeep_script.write("2> report.log \\\n")
    mirdeep_script.close()
    os.chdir("..")

