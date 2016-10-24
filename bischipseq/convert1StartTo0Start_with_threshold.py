#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 29/05/2014

import argparse

parser = argparse.ArgumentParser(description='Converts one-start bedgraph (or BED) files to zero-start bedgraph (or BED) files, with threshold for number of counts.')
parser.add_argument("-o", "--one_start_file", help="One start bedgraph (or BED) file.", required=True)
parser.add_argument("-z", "--zero_start_file_with_threshold", help="Zero start bedgraph (or BED) file.", required=True)
parser.add_argument("-t", "--threshold", help="Threshold for number of counts (methylated and unmethylated). DEFAULT=5", default=5)
args = parser.parse_args()

threshold = int(args.threshold)

one_start_file = open(args.one_start_file)
zero_start_file_with_threshold = open(args.zero_start_file_with_threshold, "w")

for one_start_line in one_start_file:
    fields = one_start_line.split()
    unmethylated_counts = int(fields[4])
    methylated_counts = int(fields[5])
    counts = unmethylated_counts + methylated_counts
    if counts >= threshold:
        chromosome = fields[0]
        start = str((int)(fields[1]) - 1) #Subtract 1 from start to convert to 0 based start format.
        zero_start_file_with_threshold.write(chromosome + "\t" + start + "\t" + "\t".join(fields[2:]) + "\n")
one_start_file.close()
zero_start_file_with_threshold.close()
