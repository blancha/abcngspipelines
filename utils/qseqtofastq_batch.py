#!/usr/bin/env python3

import os
import util

os.chdir("/sb/project/afb-431/HAI/Gen01/HAI/DNA11161/rnaseq/Libraries")

# Create script files.
for directory in os.listdir():
    scriptName = 'qseqtofastq_' + directory.replace("-", "_") + '.sh'
    script = open(scriptName, 'w')
    
    script.write("#PBS -j oe\n")
    script.write("#PBS -V\n")
    script.write("#PBS -A afb-431-ac\n\n")

    script.write("cd $PBS_O_WORKDIR\n\n")

    script.write("qseqtofastq.py " + "\\\n")
    script.write("--directory " + directory + " \\\n")
    script.write("&> " + scriptName + ".log")    
