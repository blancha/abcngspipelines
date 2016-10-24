#!/usr/bin/env python3

# Author: Alexis Blanchet-Cohen
# Date: 24/02/2014

"""Submits all scripts with .sh extension in current directory to queue."""

import argparse
import glob
import os
import os.path
import pandas
import subprocess
import socket
import sys

# Parse arguments
parser = argparse.ArgumentParser(description="Submits scripts specified in command line to queue. If no scripts are specified, submits all scripts with .sh extension in current directory to queue.")
parser.add_argument("scripts", nargs='*', help='Scripts to submit to queue DEFAULT=all files ending in .sh')
parser.add_argument("-i", "--inputDirectory", help="Input directory with BASH script files. DEFAULT=\".\"", default=".")
parser.add_argument("-d", "--dependencies", help="Hold jobs until specified job IDs have completed succesfully. E.g. --dependencies=17308513:17308515 or --dependencies=jobids.txt") 
args = parser.parse_args()

inputDirectory = os.path.abspath(args.inputDirectory)
os.chdir(inputDirectory)

if len(args.scripts) != 0:
    scripts = args.scripts
else:
    scripts = glob.glob('*.sh')

# Exit program if there are no scripts to run.
if len(scripts) == 0:
    print("There are no scripts to run.")
    sys.exit(0)

# Dependencies treatment and job submission on Guillimin.
# Create a log the IDs of the jobs submitted (not on Gen01 where screen does not return a job id)
if (socket.gethostname() != "gen01.ircm.priv"):
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

# Job submission on Gen01
else:
    # Submit jobs to queue
    for script in scripts:
        # Command to submit jobs specific to Gen01 (screen)
        command = "chmod u+x " + script + ";screen -Sdm " + script + " ./" + script # Gen01
        print(command.split(sep=";")[1])
        jobID = subprocess.check_output(command, shell=True, universal_newlines=True)
