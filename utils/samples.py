#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 24/02/2014

import argparse
import glob
import operator
import os
import pandas
import util

parser = argparse.ArgumentParser(description='Generate samples.txt file.')
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files. DEFAULT=../data/FASTQ_files/untrimmed", default="../data/FASTQ_files/untrimmed")
parser.add_argument("-l", "--lanes", help="lanes in second field of FASTQ names. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no")
parser.add_argument("-s", "--single", help="Single end sequencing. DEFAULT=no", choices=["yes", "no", "y", "n"], default="no")
parser.add_argument("-r", "--replicates", help="The group will be determining by removing the last character from the replicate name.", choices=["yes", "no", "y", "n"], default="no")
args = parser.parse_args()

# Process command-line arguments
inputDirectory = os.path.abspath(args.inputDirectory)

if (args.lanes.lower() == "yes") | (args.lanes.lower() == "y"):
    lanes = True
else:
    lanes = False

if (args.single.lower() == "yes") | (args.single.lower() == "y"):
    singleEnd = True
else:
    singleEnd = False

# Store the list of files with the extensions fastq or fastq.gz 
files = glob.glob(inputDirectory + "/*.fastq") + glob.glob(inputDirectory + "/*.fastq.gz")
files.sort()

if len(files) == 0:
    print("There are no files with the extension fastq or fastq.gz in the input directory.")
    print("Input directory: " + inputDirectory)
    exit(1)

samples = pandas.DataFrame()

if not singleEnd:
    # Read the files 2 by 2 for paired end sequencing
    for file1,file2 in util.pairwise(files):
        file1_basename = os.path.basename(file1)
        file2_basename =os.path.basename(file2)
        fields = file1_basename.split(".")
        sampleName = ""
        if os.path.splitext(file1_basename)[1] == ".gz":
            sampleName = fields[-3].split("_")[-2]
        elif os.path.splitext(file1_basename)[1] == ".fastq":
            sampleName = fields[-2].split("_")[-2] 
        else:
            sampleName = os.path.splitext(file1_basename)[1]
        sampleName = sampleName.replace("-", "_")
        if args.replicates == "yes" or args.replicates == "y":
            group = sampleName[:-1]
        else:
            group = sampleName
        if lanes:
            lane = fields[2]
        if lanes:
            samples = samples.append({"file_r1": file1, "file_r2": file2, "lane": lane, "sample": sampleName, "group": group}, ignore_index=True)
        else:
            samples = samples.append({"file_r1": file1, "file_r2": file2, "sample": sampleName, "group": group}, ignore_index=True)

    # Sort samples dataframe by samples name
    samples = samples.sort_values(by="sample")
    # Convert samples and groups columns to lower case
    samples["sample"] = samples["sample"].str.lower()
    samples["group"] = samples.group.str.lower()

    # Write samples file 
    if lanes:
        samples.to_csv("samples.txt", sep="\t", index=False, columns=["file_r1", "file_r2", "lane", "sample", "group"])
    else:
        samples.to_csv("samples.txt", sep="\t", index=False, columns=["file_r1", "file_r2", "sample", "group"])

if singleEnd:
    # Read the files one by one for single end sequencing
    for file1 in files:
        file1=os.path.basename(file1)
        fields = file1.split(".")
        sampleName = ""
        if os.path.splitext(file1)[1] == ".gz":
            sampleName = fields[-3].split("_")[-2]
        elif os.path.splitext(file1)[1] == ".fastq":
            sampleName = fields[-2].split("_")[-2]
        else:
            sampleName = os.path.splitext(file1)[1]
        sampleName = sampleName.replace("-", "_")
        if args.replicates == "yes" or args.replicates == "y":
            group = sampleName[:-1]
        else:
            group = sampleName
        if lanes:
            lane = fields[2]
        if lanes:
            samples = samples.append({"file_r1": file1, "lane": lane, "sample": sampleName, "group": group}, ignore_index=True)
        else:
            samples = samples.append({"file_r1": file1, "sample": sampleName, "group": group}, ignore_index=True)

    # Sort samples dataframe by samples name
    samples = samples.sort("sample")


    # Write samples file 
    if lanes:
        samples.to_csv("samples.txt", sep="\t", index=False, columns=["file_r1", "lane", "sample", "group"])
    else:
        samples.to_csv("samples.txt", sep="\t", index=False, columns=["file_r1", "sample", "group"])
