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


#Read the one-line TSV summary
with open(infile, "r") as f:
    line = f.readline().strip()

#TSV format: sample_id  raw_count  filtered_count  n_contigs  total_bp  strain
#Above is file format we expect, so this is how we will split the lines
sample, raw_n, filt_n, n_contigs, total_bp, strain = line.split("\t")

#Write formatted report lines
with open(outfile, "w") as out:
    out.write(f"Sample {sample} had {raw_n} read pairs before and {filt_n} read pairs after Bowtie2 filtering.\n")
    out.write(f"In the assembly of sample {sample}, there are {n_contigs} contigs > 1000 bp and {total_bp} total bp.\n")
    out.write(f"Most likely strain for sample {sample}: {strain}\n")