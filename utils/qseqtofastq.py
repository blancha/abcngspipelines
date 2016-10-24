#!/usr/bin/env python3

import argparse
import glob
import gzip
import os
import subprocess

parser = argparse.ArgumentParser(description='Converts QSEQ format to FASTQ format.')
parser.add_argument("-d", "--directory", help="Directory.", required=True)
args = parser.parse_args()

os.chdir(args.directory)

for directory in os.listdir():
    os.chdir(directory + "/Raw")
    
    # FASTQ file 1
    fastq_file_1 = open(args.directory.replace("-", "_") + "_R1.fastq", "w")
    for file in sorted(glob.glob("????1?????????????????")):
        f = gzip.open(file)
        for qseq_line in f:            
            fields = qseq_line.decode().split()
            machine = fields[0]
            run = fields[1]
            lane = fields[2]
            tile = fields[3]
            x = fields[4]
            y = fields[5]
            index = fields[6]
            read = fields[7]
            sequence = fields[8]
            quality = fields[9]
            filter = fields[10]
            if filter != 0:
                filter = "Y"
                fastq_line_1 = "@" + machine + ":" + run + ":X:" + lane + ":" + tile + ":" + x + ":" + y + " " + "Y:0:" + index
                fastq_line_2 = sequence
                fastq_line_3 = "+"
                fastq_line_4 = quality
                fastq_file_1.write(fastq_line_1 + "\n" + fastq_line_2 + "\n" + fastq_line_3 + "\n" + fastq_line_4 + "\n")
    fastq_file_1.close()
    
    # FASTQ file 2
    fastq_file_2 = open(args.directory.replace("-", "_") + "_R2.fastq", "w")
    for file in sorted(glob.glob("????2?????????????????")):
        f = gzip.open(file)
        for qseq_line in f:            
            fields = qseq_line.decode().split()
            machine = fields[0]
            run = fields[1]
            lane = fields[2]
            tile = fields[3]
            x = fields[4]
            y = fields[5]
            index = fields[6]
            read = fields[7]
            sequence = fields[8]
            quality = fields[9]
            filter = fields[10]
            if filter != 0:
                filter = "Y"
                fastq_line_1 = "@" + machine + ":" + run + ":X:" + lane + ":" + tile + ":" + x + ":" + y + " " + "Y:0:" + index
                fastq_line_2 = sequence
                fastq_line_3 = "+"
                fastq_line_4 = quality
                fastq_file_2.write(fastq_line_1 + "\n" + fastq_line_2 + "\n" + fastq_line_3 + "\n" + fastq_line_4 + "\n")
    fastq_file_2.close()