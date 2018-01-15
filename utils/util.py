#!/usr/bin/env python3

from collections import OrderedDict
from configparser import ConfigParser, ExtendedInterpolation
import os
import pandas
import socket
import subprocess
import sys

# Write header
def writeHeader(script, config, section):
    if config.has_section(section):
        if config.has_option(section, "ppn"):
            script.write("#PBS -l nodes=1:ppn=" + config.get(section, "ppn") + "\n")
        if "hydra" in socket.gethostname():
            if config.has_option(section, "mem"):
                script.write("#PBS -l mem=" + config.get(section, "mem") + "\n")
            if config.has_option(section, "vmem"):
                script.write("#PBS -l vmem=" + config.get(section, "vmem") + "\n")
        if config.has_option(section, "walltime"):
            script.write("#PBS -l walltime=" + config.get(section, "walltime") + "\n")
        if config.has_option(section, "queue"):
            script.write("#PBS -q " + config.get(section, "queue") + "\n")
        if config.has_option(section, "processors"):
            processors = config.get(section, 'processors')
        if "hydra" in socket.gethostname():
            script.write("#PBS -l epilogue=/mnt/KLEINMAN_BACKUP/home/alexis.blanchetcohen/epilogue.sh\n")
    script.write("#PBS -j oe\n")
    if not "hydra" in socket.gethostname():
        if (config.has_option('user', 'RAPid')) and (config.get('user', 'RAPid') != ""):
            script.write("#PBS -A " + config.get('user', 'RAPid') + "\n")
        else:
            script.write("#PBS -A "  + "\n")
            print("Script: " + script.name + "\nWarning: No RAPid specified.\nPlease add your RAPid to the user section in the configuration file and rerun this script.\nOr, add the RAPid manually to each script submitted to the queue.\nE.g.:\n[user]\nRAPid=xyz-123-abc", file=sys.stderr) 
    if config.has_option('user', 'email') and config.get('user', 'email') and config.has_option('user', 'send_email') and config.get('user', 'send_email') != "":
        script.write("#PBS -m " + config.get('user', 'send_email') + "\n")
        script.write("#PBS -M " + config.get('user', 'email') + "\n")
    script.write("cd $PBS_O_WORKDIR\n")
    if config.getboolean("server", "purge_modules"):
        script.write("module purge" + "\n")
    script.write("\n")

def readConfigurationFiles():
    config=ConfigParser(os.environ, interpolation=ExtendedInterpolation())
    pipelines_folder = os.environ["PIPELINES_FOLDER"]
    configuration_genomes = os.path.join(pipelines_folder, "configuration_files", "configuration_genomes.txt")
    configuration_server = os.path.join(pipelines_folder, "configuration_files", "configuration_server.txt")
    configuration_programs = os.path.join(pipelines_folder, "configuration_files", "configuration_programs.txt")
    configuration_project = os.path.join(pipelines_folder, "configuration_files", "configuration_project.txt")
    config.read([configuration_server, configuration_programs, configuration_genomes,  configuration_project, os.path.expanduser('~/configuration_user.txt'), 'configuration_project.txt'])    
    return config

def pairwise(iterable):
    "s -> (s0,s1), (s2,s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)

def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return loc

# cd to the main scripts directory, if it exists and can be located.
def cdMainScriptsDirectory():
    # Check if this is the scripts directory.
    if os.path.basename(os.getcwd()) != "scripts":
        # If it isn't check if the parent directory is the scripts directory
        if os.path.basename(os.path.dirname(os.getcwd())) == "scripts":
            # Switch to the parent directory if it is the scripts directory
            os.chdir("..")  

# Check if the inputDirectory exists, and is a directory.
def checkInputDirectory(inputDirectory):
    if not os.path.exists(inputDirectory):
        if inputDirectory == "../results/tophat":
            sys.stderr.write("Error!\n")
            sys.stderr.write("The input directory (default) does not exist.\n")
            sys.stderr.write("inputDirectory (default): " + inputDirectory + "\n")
        else:
            sys.stderr.write("Error!\n")
            sys.stderr.write("The input directory does not exist.\n")
            sys.stderr.write("inputDirectory: " + inputDirectory + "\n")
        sys.exit(1)

# Create script directory, if it does not exist yet, and cd to it.
def makeDirectory(directory, recursive=False):
    if os.path.exists(directory):
        print(directory + " already exists. No need to create.")
    else:
        if recursive:
            os.makedirs(directory)
        else:
            os.mkdir(directory)
        print("Created " + directory + ".")

# Read samples file. Create file first if it doesn't exist.
def readSamplesFile(inputDirectory=None):
    if not os.path.exists("samples.txt"):
        if inputDirectory:
            subprocess.call("samples.py --inputDirectory " + inputDirectory, shell=True)
        else:
            subprocess.call("samples.py", shell=True)
    samplesFile = pandas.read_csv("samples.txt", delim_whitespace=True)
    # Convert column names to lower case.
    samplesFile.columns = map(str.lower, samplesFile.columns)
    # Convert lanes column to string if if exists.
    if "Lane" in samplesFile.columns:
       samplesFile["Lane"] = samplesFile["Lane"].apply(str) 
    return samplesFile

# Return samples.
def getSamples(lanes=False):
    samplesFile = readSamplesFile()
    if lanes:
        samples = samplesFile["sample"] + "_lane_" + samplesFile["lane"] 
    else:
        samples = samplesFile["sample"]
        # Remove duplicate samples (Order is lost.)
        samples = list(set(samples))
    return samples
