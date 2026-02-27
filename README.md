# COMP383 Pipeline Project

This repository contains a Snakemake workflow that:
1) Downloads paired-end reads for 4 SRA samples
2) Maps reads to an HCMV reference with Bowtie2
3) Keeps only mapped reads and assembles them with SPAdes
4) Computes basic assembly statistics (contigs >1000 bp, total bp)
5) BLASTs the longest contig against a Betaherpesvirinae BLAST database to infer the most likely strain
6) Writes all results to `results/AbdulMuiz_PipelineReport.txt`

## Samples
SRR5660030, SRR5660033, SRR5660044, SRR5660045

## Software requirements
This pipeline expects the following tools to be installed and available in PATH:

- snakemake
- fasterq-dump (SRA Toolkit)
- bowtie2
- samtools
- spades.py
- blastn, makeblastdb (BLAST+)
- NCBI datasets CLI (`datasets`)
- python3 with Biopython installed

## Directory layout
- `data/raw/` : input FASTQ.gz files (downloaded from SRA)
- `ref/` : reference genomes and BLAST databases
- `results/` : all outputs (ignored by git)
- `scripts/` : helper Python scripts used by rules
- `Snakefile` : workflow definition

## Setup (download inputs)

### 1) HCMV reference genome
Place the HCMV reference FASTA here:
- `ref/hcmv.fna`


### 2) Download SRA reads (paired-end)

```bash
mkdir -p data/raw

fasterq-dump SRR5660030 -O data/raw --split-files
fasterq-dump SRR5660033 -O data/raw --split-files
fasterq-dump SRR5660044 -O data/raw --split-files
fasterq-dump SRR5660045 -O data/raw --split-files

gzip -f data/raw/*.fastq
```

### 3) Run the pipeline

```bash
snakemake --cores 4
```