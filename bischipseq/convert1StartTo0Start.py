#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 29/05/2014

import argparse
import os
import subprocess
import util

config = util.readConfigurationFiles()

parser = argparse.ArgumentParser(description='Converts one-start bedgraph (or BED) files to zero-start bedgraph (or BED) files.')
parser.add_argument("-o", "--one_start_file", help="One start bedgraph (or BED) file.", required=True)
parser.add_argument("-z", "--zero_start_file", help="Zero start bedgraph (or BED) file.", required=True)
args = parser.parse_args()

one_start_file = open(args.one_start_file)
zero_start_file = open(args.zero_start_file, "w")

for one_start_line in one_start_file:
    fields = one_start_line.split()
    chromosome = fields[0]
    start = str((int)(fields[1]) - 1) #Subtract 1 from start to convert to 0 based start format.
    zero_start_file.write(chromosome + "\t" + start + "\t" + "\t".join(fields[2:]) + "\n")

one_start_file.close()
zero_start_file.close()
