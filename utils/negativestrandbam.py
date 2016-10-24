#!/usr/bin/env python3

# Version 1.1
# Author Alexis Blanchet-Cohen
# Date: 13/04/2014

import argparse
import configparser
import glob
import os
import os.path
import shutil
import subprocess
import sys
import tempfile
import util

parser = argparse.ArgumentParser(description='Splits BAM files by strand.')
parser.add_argument("-i", "--input_bam", help="Bam file to split.")
parser.add_argument("-n", "--negative_bam", help="BAM file with reads from negative strand.")
args = parser.parse_args()

# Process command line arguments.
input_bam = os.path.abspath(args.input_bam)
negative_bam = os.path.abspath(args.negative_bam)

# Check if the input file exists.
if not os.path.exists(input_bam):
    print("Input BAM file does not exist.")
    exit(1)

# Negative BAM
temp_directory = tempfile.mkdtemp(dir=".")
header_file = os.path.join(temp_directory, "negative_header.sam")

header_command = "samtools view -H " + input_bam + " > " + header_file
negative_strand_command = "samtools view test.bam | grep 'XS:A:-' | cat header.sam - | samtools view -bS - > test2.bam"
merge_neg_strand_command = "samtools merge -f " + negative_bam + " " + neg96_bam + " " + neg144_bam
index_command = "samtools index " + negative_bam
subprocess.call (header_command, shell=True)
subprocess.call (neg144_strand_command, shell=True)
subprocess.call (merge_neg_strand_command, shell=True)
subprocess.call(index_command, shell=True)

shutil.rmtree(temp_directory)

