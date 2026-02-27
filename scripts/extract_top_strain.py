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

#Read first BLAST hit line (best hit)
with open(infile, "r") as f:
    first_line = f.readline().strip()

#Split TSV columns
fields = first_line.split("\t")

#stitle is column 3 (index 2)
stitle = fields[2]

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