#!/usr/bin/env python3
# Author Alexis Blanchet-Cohen

import argparse
import os
import os.path
import subprocess

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates the pipeline")
parser.add_argument("-t", "--type", help="Type of analysis", default="rnaseq", choices=["rnaseq", "chipseq"])
args = parser.parse_args()

type = args.type

if not os.path.basename(os.getcwd()) == "scripts":
    print("This script should only be called from the scripts folder.")
    exit(1)

if type == "rnaseq":
    print("Running tophat.py ...")
    subprocess.call("tophat.py", shell=True)
    print("Running featurecounts.py ...")
    subprocess.call("featurecounts.py", shell=True)
    print("Running cuffdiff.py ...")
    subprocess.call("cuffdiff.py", shell=True)
    print("Running bamtobigwig.py ...")
    subprocess.call("bamtobigwig.py", shell=True)

    pipeline = open("rnaseqpipeline.txt", "w")
    pipeline.write("tophat\n")
    pipeline.write("featurecounts dependencies=tophat\n")
    pipeline.write("cuffdiff dependencies=tophat\n")
    pipeline.write("cummerbund dependencies=cuffdiff\n")
