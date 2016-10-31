#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date last modified: 31/10/2016

import argparse
import fileinput
import os
import os.path
import pandas
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Downloads FASTQ files from Nanuq.")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="yes")
args = parser.parse_args()

# Read configuration files
config = util.readConfigurationFiles()
J_USER = config.get("nanuq", "J_USER")
J_PASS = config.get("nanuq", "J_PASS")

command = "unzip -o readSetLinks.zip"
subprocess.call(command, shell=True)

f_open = open("run_wget.sh", "r")
lines = f_open.readlines()
lines[2] = lines[2].replace("J_USER=\"\"","J_USER=" + J_USER)
lines[3] = lines[3].replace("J_PASS=\"\"","J_PASS=" + J_PASS)
del(lines[5:8])
lines.insert(1, "#PBS -A feb-684-ac\ncd $PBS_O_WORKDIR\n\n")
f_open.close()

f_open = open("run_wget.sh", "w")
f_open.write(''.join(lines))
f_open.close()

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py run_wget.sh", shell=True)
