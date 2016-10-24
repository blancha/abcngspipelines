#!/usr/bin/env python3
# Author Alexis Blanchet-Cohen

import argparse
import os.path
import subprocess

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Submits the jobs found in the pipeline file to the queue.")
parser.add_argument("pipeline", nargs=1, help="Pipeline file. DEFAULT=./pipeline.txt", default="pipeline.txt")
args = parser.parse_args()

pipelineFile = args.pipeline[0]

if not os.path.exists(pipelineFile):
    print("Specified pipeline file does not exist:")
    print(pipelineFile)
    exit(1)

# Read pipelines.txt line by line
for line in open(pipelineFile):
    line = line.lower()
    fields=line.split("dependencies=")
    jobs_folder=fields[0].strip()
    command = "submitJobs.py " + jobs_folder
    if len(fields) > 1:
        dependencies=fields[1].strip()
        command += " --dependencies=" + dependencies
        #subprocess.call(command, shell=True)
    print(command)
    subprocess.call(command, shell=True)
