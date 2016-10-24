#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 24/02/2014

import argparse
import os
import shutil

parser = argparse.ArgumentParser(description='Make project directory structure.')
parser.add_argument("-n", "--name", help="Name of the project.", required=True)
parser.add_argument("-c", "--config_file", help="Configuration file to be copied to Analysis/scripts/config_project.txt")
args = parser.parse_args()

# Get absolute path to configuration file.
if args.config_file != None:
    config_file = os.path.abspath(args.config_file)

# Create project directories, if they do not exist yet, and change to main project directory.
if not os.path.exists(args.name):
    os.mkdir(args.name)
os.chdir(args.name)

# Create output directories, if they do not exist yet..
if not os.path.exists("FASTQ/untrimmed"):
    os.makedirs("FASTQ/untrimmed")

if not os.path.exists("Analysis/scripts"):
    os.makedirs("Analysis/scripts")

shutil.copy(config_file, "Analysis/scripts/config_project.txt")

