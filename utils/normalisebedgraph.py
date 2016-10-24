#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 12/04/2014

import argparse
import os
import subprocess
import util

parser = argparse.ArgumentParser(description='Normalises bedgraph files.')
parser.add_argument("-b", "--before_normalisation", help="Bedgraph file before normalisation.", required=True)
parser.add_argument("-a", "--after_normalisation", help="Bedgraph file after normalisaion.", required=True)
parser.add_argument("-s", "--scaling_factor", help="Scaling factor.", required=True)
args = parser.parse_args()

if args.scaling_factor == 1:
    subprocess.call("cp "+ args.before_normalisation + " " + args_after_normalisation, shell=True)

bedgraph_before_normalisation = open(args.before_normalisation)
bedgraph_after_normalisation = open(args.after_normalisation, "w")

for before_normalisation_line in bedgraph_before_normalisation:
    fields = before_normalisation_line.split()
    normalised_value = "{0:.2f}".format(float(args.scaling_factor) * float(fields[3]))
    after_normalisation_line = fields[0] + "\t" + fields[1] + "\t" + fields[2] + "\t" + normalised_value + "\n"
    bedgraph_after_normalisation.write(after_normalisation_line)

bedgraph_before_normalisation.close()
bedgraph_after_normalisation.close()
