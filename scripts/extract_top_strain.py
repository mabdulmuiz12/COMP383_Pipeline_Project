#!/usr/bin/env python3

import sys
import argparse


#function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="<add description of what script does>")
    parser.add_argument("-i", "--input", 
    help="input file", 
    required=True)
    parser.add_argument("-o", "--output",
    help="output file",
    required=True)
    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
infile = arguments.input
outfile = arguments.output

#Read the first (top) BLAST hit line
with open(infile, "r") as f:
    line = f.readline().strip()

#Split the tab-delimited BLAST outfmt 6 line into columns
cols = line.split("\t")

#stitle is the LAST column in our BLAST outfmt
stitle = cols[-1]

#Default name in case parsing fails
name = "UNKNOWN"

#Try extracting after "strain "
key1 = "strain "
idx1 = stitle.lower().find(key1)
if idx1 != -1:
    after = stitle[idx1 + len(key1):]
    name = after.split(",", 1)[0].strip()
else:
    #If no "strain", then try extracting after "clone "
    key2 = "clone "
    idx2 = stitle.lower().find(key2)
    if idx2 != -1:
        after = stitle[idx2 + len(key2):]
        name = after.split(",", 1)[0].strip()

#Write extracted name
with open(outfile, "w") as out:
    out.write(f"{name}\n")