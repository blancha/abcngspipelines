#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 26/05/2014

import argparse
import os
import subprocess
import util

config = util.readConfigurationFiles()

parser = argparse.ArgumentParser(description='Extracts methylation counts from Bismark cov file, and create 2 bedgraph files, with methylated and unmethylated counts, respectively.')
parser.add_argument("-c", "--coverage", help="Bismark coverage file.", required=True)
parser.add_argument("-m", "--methylated", help="Methylated counts bedgraph file.", required=True)
parser.add_argument("-u", "--unmethylated", help="Unmethylated counts bedgraph file.", required=True)
args = parser.parse_args()

coverage = open(args.coverage)
methylated_bedgraph = open(args.methylated, "w")
unmethylated_bedgraph = open(args.unmethylated, "w")

for coverage_line in coverage:
    fields = coverage_line.split()
    chromosome = fields[0]
    start = fields[1]
    end = fields[2]
    methylated_count = fields[4]
    unmethylated_count = fields[5]
    methylated_bedgraph.write(chromosome +  "\t" + start + "\t" + end + "\t" + methylated_count + "\n")
    unmethylated_bedgraph.write(chromosome +  "\t" + start + "\t" + end + "\t" + unmethylated_count + "\n")

coverage.close()
methylated_bedgraph.close()
unmethylated_bedgraph.close()
