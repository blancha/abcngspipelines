#!/usr/bin/env python

import glob
import os
import sys
sys.path.append('/sb/software/areas/ircm/tools/ngs_pipelines/utils/')
import util

config=util.readConfigurationFiles()

# Create scripts directory, if it does not exist yet, and cd to it.
if not os.path.exists("samtoolsIndex"):
    os.mkdir("samtoolsIndex")
os.chdir("samtoolsIndex")

samplesFile = open("../samples.txt")
samplesLine = samplesFile.readlines()[1::2]

samples = []
lanes = []

for line in samplesLine:
    samples.append(line.split()[3].split("_")[-2])
    lanes.append(line.split()[2].split("_")[-1])
    
uniquesamples = list(set(samples))

for uniquesample in uniquesamples:
    # Create script
    script = open("samtoolsIndex_" + uniquesample + ".sh", "w")
    util.write_header(script, config, "samtoolsIndex")
    script.write("samtools sort ../../bismark/" + uniquesample + "/" + uniquesample + "_bismark_pe.bam ../../bismark/" + uniquesample + "/" + uniquesample + "_bismark_pe_sorted")     
    script.write ("\n\n")
    script.write("samtools index ../../bismark/" + uniquesample + "/" + uniquesample + "_bismark_pe_sorted.bam \n")    

script.close()

