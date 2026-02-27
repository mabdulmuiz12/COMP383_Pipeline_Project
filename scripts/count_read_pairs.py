import sys
import argparse
import gzip
from Bio import SeqIO

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

#Open gzipped FASTQ and count reads using Biopython
with gzip.open(infile, "rt") as handle:
    count = sum(1 for _ in SeqIO.parse(handle, "fastq"))

#Write count to output file
with open(outfile, "w") as out:
    out.write(f"{count}\n")

