#!/usr/bin/env python3

from collections import OrderedDict
import argparse
import glob
import os
import os.path
import pandas

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates file with the ucsc tracks.")
parser.add_argument("-s", "--samples", help="Samples file, with description of samples. DEFAULT=../../../scripts/samples.txt", default="../../../scripts/samples.txt")
parser.add_argument("-o", "--output", help="Output file, with ucsc tracks. DEFAULT=ucsc_tracks.txt", default="ucsc_tracks.txt")
parser.add_argument("-c", "--conditions", help="Conditions files? Generated from means of replicate samples. Not available for unnormalized bigWig files. DEFAULT=yes", choices=["yes", "no", "y", "n"], default="yes")
parser.add_argument("-r", "--stranded", help="Stranded bigWig files? DEFAULT=yes", choices=["yes", "no", "y", "n"], default="yes")
args = parser.parse_args()

# Process the command line arguments.
samplesFile = os.path.abspath(args.samples) 
outputFile = os.path.abspath(args.output)
if (args.conditions.lower() == "yes") | (args.conditions.lower() == "y"):
    conditions = True
else:
    conditions = False
if (args.stranded.lower() == "yes") | (args.stranded.lower() == "y"):
    stranded = True
else:
    stranded = False

ucsc_tracks = open(outputFile, "w")

projectDirectory = os.getcwd()
fields = projectDirectory.split("/")
lab = fields[-5]
project = fields[-4]
cwd = fields[-1]

# Read samples file.
samplesFile = pandas.read_csv(samplesFile, delim_whitespace=True)

# Strands
if stranded:
    strands = ["", "_positive", "_negative"]
else:
    strands = [""]

# Iterate over samples file
for index, row in samplesFile.iterrows():
    if conditions:
        # Test if condition is new
        if (index==0):
            new_condition = True
        else:
            if (row["condition"] != previous_condition):
                new_condition = True
            else:
                new_condition = False
        # If condition is new, generate tracks for new condition.      
        if new_condition: 
            for strand in strands:
                ucsc_tracks.write("track type=bigWig name=" + row["condition"] + strand + " ")
                if stranded & (strand != ""): 
                    ucsc_tracks.write("visibility=0")
                else:
                    ucsc_tracks.write("visibility=2")
                ucsc_tracks.write("color=0,0,0 ")
                ucsc_tracks.write("bigDataUrl=" + os.path.join("http://biosrv02.ircm.qc.ca/", lab, project, cwd, row["condition"] + strand + ".bw"))
                ucsc_tracks.write("\n")
    # Generate tracks for sample.
    for strand in strands:
        ucsc_tracks.write("track type=bigWig name=" + row["sample"] + strand + " ")
        if (strand != "") :
            ucsc_tracks.write("visibility=0 ")
        elif conditions:
            ucsc_tracks.write("visibility=0 ")
        else:
            ucsc_tracks.write("visibility=2 ")
        ucsc_tracks.write("color=0,0,0 ")
        ucsc_tracks.write("bigDataUrl=" + os.path.join("http://biosrv02.ircm.qc.ca/", lab, project, cwd, row["sample"] + strand + ".bw"))
        ucsc_tracks.write("\n")
    # Reset condition
    previous_condition = row["condition"]

ucsc_tracks.close()
                       
