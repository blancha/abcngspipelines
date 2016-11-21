#!/usr/bin/env python3

# Author: Alexis Blanchet-Cohen
# Date: 24/02/2014

"""Submits all scripts with .sh extension in current directory to queue."""

import argparse
import glob
import os.path
import pandas
import subprocess
import sys

# Parse arguments
parser = argparse.ArgumentParser(description="Submits scripts specified in command line to queue. If no scripts are specified, submits all scripts with .sh extension in current directory to queue.")
parser.add_argument("scripts", nargs='*', help='Scripts to submit to queue, or directory containing scripts DEFAULT=all files ending in .sh')
parser.add_argument("-d", "--dependencies", help="Hold jobs until specified job IDs have completed succesfully. E.g. --dependencies=17308513:17308515 or --dependencies=jobids.txt") 
args = parser.parse_args()

# Dependencies treatment
# Create a log the IDs of the jobs submitted (not on Gen01 where screen does not return a job id)
# Process dependencies if any given in argument.
if args.dependencies != None:
    if os.path.isdir(args.dependencies):
        dependenciesFile = pandas.read_csv(os.path.join(args.dependencies, "jobids.txt"), sep="\t", dtype=str)
        dependencies = dependenciesFile["jobid"]
        dependencies=":".join(dependencies)
    elif os.path.isfile(os.path.basename(args.dependencies)):
        dependenciesFile = pandas.read_csv(args.dependencies, sep="\t", dtype=str)
        dependencies = dependenciesFile["jobid"]
        dependencies=":".join(dependencies)
    else:
        dependencies=args.dependencies

# If a directory is given in argument, find all scripts in directory.
if (len(args.scripts) != 0 and os.path.isdir(args.scripts[0])):
    inputDirectory = args.scripts[0]
    os.chdir(inputDirectory)
    scripts = glob.glob("*.sh")
# If no scripts or directory are given in argument, find all scripts in current directory.        
elif (len(args.scripts) == 0):
    scripts = glob.glob("*.sh")
# Otherwise, the scripts were given in argument.
else:
    scripts = args.scripts

# Exit program if there are no scripts to run.
if len(scripts) == 0:
    print("There are no scripts to run.")
    sys.exit(1)

# Create job id log file
jobIDFile = open("jobids.txt", "w")
jobIDFile.write("script\tjobid\n")

# Submit jobs to queue
for script in scripts:
    command = "qsub ./" + script # Guillimin
    # Add dependencies, if any
    if (args.dependencies != None): 
        command=command + " -W depend=afterok:" + dependencies
    print(command)
    jobID = subprocess.check_output(command, shell=True, universal_newlines=True)
    # Write the job IDS to a log file (not on Gen01 where screen does not return a job id)
    print(jobID, sep="")
    jobIDFile.write(script + "\t" + jobID.split(sep=".")[0] + "\n")

