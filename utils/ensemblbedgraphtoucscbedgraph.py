#!/usr/bin/env python3

# Version 1.0
# Author Alexis Blanchet-Cohen
# Date: 12/04/2014

import argparse
import os
import pandas
import util

config = util.readConfigurationFiles()

parser = argparse.ArgumentParser(description='Converts bedgraph files in Ensembl format to UCSC format.')
parser.add_argument("-e", "--ensembl_bedgraph", help="Ensembl bedgraph.", required=True)
parser.add_argument("-u", "--ucsc_bedgraph", help="Generated UCSC bedgraph.", required=True)
parser.add_argument("-d", "--dictionary", help="Tab-delimited text file with the Ensembl symbols and the corresponding UCSC symbols. Default=/sb/project/afb-431-ac/genomes/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/ucsctoEnsembl.txt", default="/sb/project/afb-431-ac/genomes/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/ucsctoEnsembl.txt")
args = parser.parse_args()

ensembl_bedgraph = open(args.ensembl_bedgraph)
ucsc_bedgraph = open(args.ucsc_bedgraph, "w")
ucsctoEnsembl_file = pandas.read_csv(args.dictionary, sep="\t")

# Build dictionary
ucsctoEnsembl_dictionary = ucsctoEnsembl_file.set_index("ensembl")["ucsc"].to_dict()

for ensembl_line in ensembl_bedgraph:
    fields = ensembl_line.split()
    ensembl_chromosome = fields[0]
    if ensembl_chromosome in ucsctoEnsembl_dictionary:
        ucsc_line = ucsctoEnsembl_dictionary[ensembl_chromosome] + "\t" + fields[1] + "\t" + fields[2] + "\t" + fields[3] + "\n"
        ucsc_bedgraph.write(ucsc_line)

ensembl_bedgraph.close()
ucsc_bedgraph.close()
