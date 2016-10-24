#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 04/06/2014

import argparse

parser = argparse.ArgumentParser(description='Calculates average methylation level across regions.')
parser.add_argument("-i", "--input_file", help="Input BED file.", required=True)
parser.add_argument("-o", "--output_file", help="Output text file with average methylation level across all regions.", required=False, default="results.txt")

args = parser.parse_args()

input_file = open(args.input_file)
output_file = open(args.output_file, "w")

total_bases_covered = 0
total_average_methylation_percentages = 0

for input_line in input_file:
    fields = input_line.split()
    start = int(fields[1])
    end = int(fields[2])
    methylation_percentage = float(fields[4])
    bases_covered = end - start
    average_methylation_percentage = methylation_percentage * bases_covered 
    total_average_methylation_percentages += average_methylation_percentage
    total_bases_covered += bases_covered 

average_methylation_percentage = total_average_methylation_percentages / total_bases_covered

output_file.write("Average methylation (%)\n")
output_file.write("{0:.2f}".format(average_methylation_percentage))

input_file.close()
output_file.close()
