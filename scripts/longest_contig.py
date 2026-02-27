import sys
import argparse
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

#store the longest contig record found so far
longest_record = None 

#loop through contigs and keep the longest one
for record in SeqIO.parse(infile, "fasta"):
    if (longest_record is None) or (len(record.seq) > len(longest_record.seq)):
        longest_record = record

#write the longest contig to a FASTA output file
with open(outfile, "w") as out:
    SeqIO.write(longest_record, out, "fasta")
