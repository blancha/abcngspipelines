#!/bin/bash

date > softwareVersions.txt

echo -e "\n" >> softwareVersions.txt
echo "TopHat2 " >> softwareVersions.txt
which tophat2 >> softwareVersions.txt

echo -e "\n" >> softwareVersions.txt
echo "Bowtie2 " >> softwareVersions.txt
which bowtie2 >> softwareVersions.txt

echo -e "\n" >> softwareVersions.txt
echo "Cufflinks " >> softwareVersions.txt
which cufflinks >> softwareVersions.txt

echo -e "\n" >> softwareVersions.txt
echo "Java Classpath: " >> softwareVersions.txt
echo $CLASSPATH >> softwareVersions.txt

echo -e "\n" >> softwareVersions.txt
echo "Python path: " >> softwareVersions.txt
echo $PYTHONPATH >> softwareVersions.txt