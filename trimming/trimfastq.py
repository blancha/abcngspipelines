#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 09/06/2014

import argparse
import glob
import os
import subprocess
import util

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates trimFastq.py scripts.")
parser.add_argument("-s", "--scriptsDirectory", help="Scripts directory.", default="trimfastq")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files.", default="../data/FASTQ_files/untrimmed")
parser.add_argument("-o", "--outputDirectory", help="Output directory with trimmed FASTQ files.", default="../data/FASTQ_files/trimmed")
parser.add_argument("-q", "--submitJobsToQueue", help="Submit jobs to queue immediately.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# If not in the main scripts directory, cd to the main scripts directory, if it exists.
util.cdMainScriptsDirectory()

# Process the command line arguments.
scriptsDirectory = os.path.abspath(args.scriptsDirectory) 
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Check if the inputDirectory exists, and is a directory.
util.checkInputDirectory(args.inputDirectory)

# Read configuration files
config = util.readConfigurationFiles()

header = config.getboolean("server", "PBS_header")
minlength=config.get("trimfastq", "minlength")

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists(scriptsDirectory):
    os.mkdir(scriptsDirectory)
os.chdir(scriptsDirectory)

# Create output directory, if it does not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Store the list of files with the extensions fastq or fastq.gz 
files = glob.glob(inputDirectory + "/*.fastq") + glob.glob(inputDirectory + "/*.fastq.gz")

# Create script files.
for file in files:
    file = os.path.basename(file)
    scriptName = 'trimfastq_' + file + '.sh'
    script = open(scriptName, 'w')
    if header:
        util.writeHeader(script, config, "trimfastq")
    script.write("#Deactivate Python 3 virtual environment, and activate Python 2 virtual environment to be able to use TopHat.\n")
    virtual_env_directory = os.path.dirname(os.environ["VIRTUAL_ENV"])
    script.write("source " + os.path.join(virtual_env_directory, "python2.7/bin/activate") + "\n\n")    
    script.write("trimFastq.py"   + " \\\n")
    script.write(os.path.relpath(os.path.join(inputDirectory, file)) + " \\\n")
    script.write(os.path.relpath(os.path.join(outputDirectory, file)) + " \\\n")
    script.write(minlength + " \\\n")
    script.write("&> " + scriptName + ".log")
    script.close()    

if (args.submitJobsToQueue.lower() == "yes") | (args.submitJobsToQueue.lower() == "y"):
    subprocess.call("submitJobs.py", shell=True)
