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

lengths = []  # Store contig lengths

#Parse FASTA file
for record in SeqIO.parse(infile, "fasta"):
     #Get length of each contig
    lengths.append(len(record.seq)) 

#Filter for contigs strictly greater than 1000 bp
over_1000 = [L for L in lengths if L > 1000]

#Count contigs >1000 bp
n_contigs = len(over_1000)

#Calculate total combined length of contigs >1000 bp
total_bp = sum(over_1000)

#Write results to output file
with open(outfile, "w") as out:
    out.write(f"{n_contigs}\t{total_bp}\n")