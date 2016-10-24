#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 12/04/2014

import argparse
import math
import os
import util

parser = argparse.ArgumentParser(description='Log2 transform bedgraph files.')
parser.add_argument("-b", "--before_log2", help="Bedgraph file before normalisation.", required=True)
parser.add_argument("-a", "--after_log2", help="Bedgraph file after normalisaion.", required=True)
args = parser.parse_args()

bedgraph_before_log2 = open(args.before_log2)
bedgraph_after_log2 = open(args.after_log2, "w")

for before_log2_line in bedgraph_before_log2:
    fields = before_log2_line.split()
    log2_value = "{0:.2f}".format(math.log(float(fields[3]), 2))
    after_log2_line = fields[0] + "\t" + fields[1] + "\t" + fields[2] + "\t" + log2_value + "\n"
    bedgraph_after_log2.write(after_log2_line)

bedgraph_before_log2.close()
bedgraph_after_log2.close()
