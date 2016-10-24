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
parser.add_argument("-p", "--positive_bam", help="BAM file with reads from positive strand.")
parser.add_argument("-n", "--negative_bam", help="BAM file with reads from negative strand.")
args = parser.parse_args()

# Process command line arguments.
input_bam = os.path.abspath(args.input_bam)
positive_bam = os.path.abspath(args.positive_bam)
negative_bam = os.path.abspath(args.negative_bam)

# Check if the input file exists.
if not os.path.exists(input_bam):
    print("Input BAM file does not exist.")
    exit(1)

# Negative BAM
temp_directory = tempfile.mkdtemp(dir=".")
neg96_bam = os.path.join(temp_directory, "neg96.bam")
neg144_bam = os.path.join(temp_directory, "neg144.bam")

neg96_strand_command = "samtools view -h -f 96 " + input_bam  + " | samtools view -bS - > " + neg96_bam
neg144_strand_command = "samtools view -h -f 144 " + input_bam + " | samtools view -bS - > " + neg144_bam
merge_neg_strand_command = "samtools merge -f " + negative_bam + " " + neg96_bam + " " + neg144_bam
index_command = "samtools index " + negative_bam
subprocess.call (neg96_strand_command, shell=True)
subprocess.call (neg144_strand_command, shell=True)
subprocess.call (merge_neg_strand_command, shell=True)
subprocess.call(index_command, shell=True)

shutil.rmtree(temp_directory)

# Positive BAM
temp_directory = tempfile.mkdtemp(dir=".")
pos80_bam = os.path.join(temp_directory, "pos80.bam")
pos160_bam = os.path.join(temp_directory, "pos160.bam")

pos80_strand_command = "samtools view -h -f 80 " + input_bam  + " | samtools view -bS - > " + pos80_bam
pos160_strand_command = "samtools view -h -f 160 " + input_bam + " | samtools view -bS - > " + pos160_bam
merge_pos_strand_command = "samtools merge -f " + positive_bam + " " + pos80_bam + " " + pos160_bam
index_command = "samtools index " + positive_bam
subprocess.call (pos80_strand_command, shell=True)
subprocess.call (pos160_strand_command, shell=True)
subprocess.call (merge_pos_strand_command, shell=True)
subprocess.call(index_command, shell=True)

shutil.rmtree(temp_directory)
