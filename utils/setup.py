#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 24/02/2014

import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='Set up analysis directories and configuration files.')
parser.add_argument("directory", nargs=1, help="Project directory. If the directory does not exist, it will be created.")
parser.add_argument("-c", "--configuration", help="Configuration file. Choices=bisseq, chipseq, rnaseq, or user-defined file", default="rnaseq")
args = parser.parse_args()

pipelines_folder = os.environ["PIPELINES_FOLDER"]

# Process command line arguments
directory = args.directory[0]

if args.configuration == "bisseq":
    configuration_file = os.path.join(pipelines_folder, "configuration_files/bisseq/configuration_project.txt")
elif args.configuration == "chipseq":
    configuration_file = os.path.join(pipelines_folder, "configuration_files/chipseq/configuration_project.txt")
elif args.configuration == "rnaseq":
    configuration_file = os.path.join(pipelines_folder, "configuration_files/rnaseq/configuration_project.txt")
else:
    configuration_file = os.path.abspath(args.configuration)

# Create project directory, if they do not exist yet, and change to main project directory.
if not os.path.exists(directory):
    os.mkdir(directory)
    print("Project directory created: " + directory)
else:
    print("No need to create project directory. Directory already exists: " + directory)
os.chdir(directory)

#Create subdirectoroes, if they do not exist yet.
#if not os.path.exists("data/BAM_files"):
#    os.makedirs("data/BAM_files")

if not os.path.exists("data/FASTQ_files/untrimmed"):
    os.makedirs("data/FASTQ_files/untrimmed")

if not os.path.exists("scripts"):
    os.makedirs("scripts")

if not os.path.exists("results"):
    os.makedirs("results")

# Create link to report file
subprocess.call("ln -s scripts/report/report.pdf report.pdf", shell=True)

# Copy configuration file to Analysis/scripts
subprocess.call("cp " + configuration_file + " scripts", shell=True)
print("Configuration file copied to " + directory + "/scripts:\n" + configuration_file)


