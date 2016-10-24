import glob
import os
import subprocess

directories = os.listdir(".")

for directory in directories:
    os.chdir(directory)
    bamFile = str(glob.glob("*.bam"))
    oldName = bamFile[2:-17]
    newName = os.path.split(os.getcwd())[1]
    command = "rename " + oldName + " " + newName + " *"
    #print(command)
    subprocess.call(command, shell=True)
    os.chdir("..")
        
for directory in directories:
    os.chdir(directory)
    bamFile = str(glob.glob("*.fastq.gz_unmapped_reads_2.txt"))
    oldName = bamFile[2:-32]
    newName = os.path.split(os.getcwd())[1]
    command = "rename " + oldName + " " + newName + " *"
    #print(command)
    subprocess.call(command, shell=True)
    os.chdir("..")
    
for directory in directories:
    os.chdir(directory)
    command = "rename '.fastq.gz' '' *"
    #print(command)
    subprocess.call(command, shell=True)
    os.chdir("..")

