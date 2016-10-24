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

parser = argparse.ArgumentParser(description='Separates specified strand from BAM file.')
parser.add_argument("-i", "--input_bam", help="Bam file to split.")
parser.add_argument("-s", "--strand", help="Strand to separate. DEFAULT=positive", default="positive", choices=["positive", "negative"])
args = parser.parse_args()

# Process command line arguments.
input_bam = os.path.abspath(args.input_bam)
strand = args.strand
if (strand == "positive"):
    strand_sign = "+"
else:
    strand_sign = "-"

# Check if the input file exists.
if not os.path.exists(input_bam):
    print("Input BAM file does not exist.")
    exit(1)

# Header command
header_file = tempfile.NamedTemporaryFile()
header_command = "samtools view -H " + input_bam + " > " + header_file.name
subprocess.call(header_command, shell=True)

# Separate strand command
strand_command = "samtools view " + input_bam + " | grep 'XS:A:" + strand_sign + "' | cat " + header_file.name + " - | samtools view -bS - > " + input_bam.replace(".bam", "_" + args.strand + ".bam")
subprocess.call(strand_command, shell=True)

# Index command
index_command = "samtools index " + input_bam.replace(".bam", "_" + strand + ".bam")
subprocess.call(index_command, shell=True)

