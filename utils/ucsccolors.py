#!/usr/bin/env python3

import ast
import argparse

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates file with the ucsc tracks.")
parser.add_argument("-d", "--dictionary", help="Dictionary of colors. DEFAULT={'mbne': '0,0,255', 'memb': '0,0,255', 'nucl': '0,255,0', 'cyto': '255,0,0', 'total': '0,0,0'}", default="{'mbne': '0,0,255', 'memb': '0,0,255', 'nucl': '0,255,0', 'cyto': '255,0,0', 'total': '0,0,0'}")
args = parser.parse_args()

# Process the command line arguments.
dictionary = ast.literal_eval(args.dictionary)

old_ucsc_tracks = open("ucsc_tracks.txt").readlines()

new_ucsc_tracks = open("ucsc_tracks.txt", "w")

for old_line in old_ucsc_tracks:
    fields = old_line.split()
    name = fields[2].split("name=")[1].lower()
    old_color = fields[4].split("color=")[1].lower()
    for key, value in dictionary.items():   # iter on both keys and values
        if key in name:
            new_color = dictionary[key]
            break
    new_line = old_line.replace(old_color, new_color)
    new_ucsc_tracks.write(new_line)

new_ucsc_tracks.close()
