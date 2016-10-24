#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 12/04/2014

import argparse
import math
import os
import os.path
import subprocess
import tempfile
import util

parser = argparse.ArgumentParser(description="Computes means of bedgraph files.")
parser.add_argument("-b", "--bedgraphs", help="Comma-separated list of bedgraph files from which mean should be computed. Example, file1.bedgraph, file2.bedgraph", required=True, type=str)
parser.add_argument("-o", "--output", help="Bedgraph file with means of the bedgraph files given in input.", required=True)
args = parser.parse_args()

# Process command-line arguments
bedgraphs = args.bedgraphs.replace(" ", "").split(",")
bedgraphs_mean = open(args.output, "w")

# Generate union of bedgraphs given in argument
union_bedgraphs = "union_" 
for bedgraph in bedgraphs[:-1]:
    union_bedgraphs += os.path.splitext(os.path.basename(bedgraph))[0] + "_"
union_bedgraphs += os.path.splitext(os.path.basename(bedgraphs[-1]))[0] + ".bedgraph"
command = "bedtools unionbedg -i " 
for bedgraph in bedgraphs:
    command += os.path.relpath(bedgraph) + " " 
command  += "> " + os.path.relpath(union_bedgraphs)
print("Running: " + command)
subprocess.call(command, shell=True)
print("Completed")

# Generate bedgraph_means
for line in open(union_bedgraphs):
    fields = line.split()
    # Convert fields with coverage to floats
    coverages = [float(i) for i in fields[4:]]
    mean = sum(coverages)/len(coverages)
    bedgraphs_mean.write("\t".join(fields[0:3]) + "\t" + str(mean) + "\n")
bedgraphs_mean.close()
print("Generated: " + args.output)

os.remove(union_bedgraphs)
