#!/bin/bash

#################################################
# Parses Homer output file.                     #
#                                               #
# Generates 4 statistics files:                 #
# tss.stats.csv, exon.stats.csv                 #
# intron.stats.csv & tss.distance.csv.          #
# The 4 files are stored in the graphs folder.  #
#                                               #
# Author: Alexis Blancehet-Cohen                #
# Modified from script written by Maxime Caron  #
#################################################

if [ ! -n "$1" ]
then
  echo "Wrong number of arguments"
  echo "Usage: parseHomerOuput.sh homerOutputFile.csv"
  exit
fi  

homerFileName=$1

# Create graphs directory is it doesn't exist
if [ ! -d "graphs" ]; then
    `mkdir graphs`
fi

# Create graphs/statistics directory is it doesn't exist
if [ ! -d "graphs/statistics" ]; then
    `mkdir graphs/statistics`
fi

# Create graphs/pdf directory is it doesn't exist
if [ ! -d "graphs/pdf" ]; then
    `mkdir graphs/pdf`
fi

a=`cat $homerFileName | awk -F'\t' '{print $8}' | awk '$1 == "exon"' | wc -l`
b=`cat $homerFileName | awk -F'\t' '{print $8}' | awk '$1 == "intron"' | wc -l`
c=`cat $homerFileName | awk -F'\t' '$8 !~ "exon" && $8 !~ "intron"' | awk -F"\t" '$10>-2000 && $10<0' | wc -l`
d=`cat $homerFileName | awk -F'\t' '$8 !~ "exon" && $8 !~ "intron"' | awk -F"\t" '$10<=-2000 && $10>-10000' | wc -l`
e=`cat $homerFileName | awk -F'\t' '$8 !~ "exon" && $8 !~ "intron"' | awk -F"\t" '$10<=-10000 && $10>-100000' | wc -l`
f=`sed 1d $homerFileName | awk -F'\t' '$8 !~ "exon" && $8 !~ "intron"' | awk -F'\t' '$10<-100000 || $10>100000' | wc -l`
g=`sed 1d $homerFileName | wc -l`
h=`expr $g - $a - $b - $c - $d - $e - $f`
echo "exon,intron,proximal,distal,5d,gene_desert,other" > graphs/statistics/tss.stats.csv
echo "$a,$b,$c,$d,$e,$f,$h" >> graphs/statistics/tss.stats.csv
cat $homerFileName | awk -F'\t' '{print $8}' | awk '$1 == "exon"' | awk -F',' '{print $2}' | sed -e 's/)//g' | awk '{print $2}' > graphs/statistics/exon.stats.csv
cat $homerFileName | awk -F'\t' '{print $8}' | awk '$1 == "intron"' | awk -F',' '{print $2}' | sed -e 's/)//g' | awk '{print $2}' > graphs/statistics/intron.stats.csv
sed 1d $homerFileName | awk -F'\t' '{print $10}' > graphs/statistics/tss.distance.csv