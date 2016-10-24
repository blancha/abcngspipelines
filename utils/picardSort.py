#!/usr/bin/env python3

import glob
import os
import util

config=util.readConfigurationFiles()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists("picardSort"):
    os.mkdir("picardSort")
os.chdir("picardSort")

samplesFile = open("../samples.txt")
samplesLine = samplesFile.readlines()[1::2]

samples = []
lanes = []

for line in samplesLine:
    samples.append(line.split()[3].split("_")[-2])
    
uniquesamples = list(set(samples))

for uniquesample in uniquesamples:
    # Create script
    script = open("picardSort_" + uniquesample + ".sh", "w")
    util.write_header(script, config, "picardSort")
    script.write("samtools sort ../../bismark/" + uniquesample + "/" + uniquesample + ".bam ../../bismark/" + uniquesample + "/" + uniquesample + "_sorted")     
    script.write ("\n\n")
    script.write("samtools index ../../bismark/" + uniquesample + "/" + uniquesample + "_sorted.bam \n")    
    script.close()
