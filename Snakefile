# List of SRA sample IDs used in this project
SAMPLES = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]

#Final target of the entire thing
rule all:
    input:
        "results/PipelineReport.txt"

#Build Bowtie2 index for HCMV reference genome
#Gives us a genome index before reads can me mapped
rule bowtie2_index:
    #Input is HCMV genome FASTA
    input:
        fasta="ref/hcmv.fna"
    #Output is Bowtie2 index files
    output:
        idx="ref/bowtie2/hcmv.1.bt2"
    shell:
        """
        mkdir -p ref/bowtie2
        bowtie2-build {input.fasta} ref/bowtie2/hcmv
        """


rule bowtie2_map:
    input:
        #Bowtie2 requires the genome index prefix
        idx="ref/bowtie2/hcmv.1.bt2",
        #Paired-end FASTQ files for each sample
        r1="data/raw/{sample}_1.fastq.gz",
        r2="data/raw/{sample}_2.fastq.gz"
    output:
         #SAM alignment output file
        sam="results/mapped/{sample}.sam",
         #Log file capturing error and output messages from Bowtie2 into a file
        log="logs/bowtie2/{sample}.log"
    shell:
        r"""
        #Create output directories if they do not exist
        mkdir -p results/mapped logs/bowtie2

        #Run Bowtie2 alignment
        bowtie2 -x ref/bowtie2/hcmv \
            -1 {input.r1} \
            -2 {input.r2} \
            -S {output.sam} \
            2> {output.log}
        """

rule sam_to_mapped_fastq:
    input:
        sam="results/mapped/{sample}.sam"
    output:
        out_r1="results/filtered/{sample}_1.fastq.gz",
        out_r2="results/filtered/{sample}_2.fastq.gz"
    threads: 2
    shell:
        r"""
        #Making sure output directory exists
        mkdir -p results/filtered

        # Convert SAM to BAM and remove unmapped reads
        #And convert remaining mapped reads back to paired FASTQ
        samtools view -@ {threads} -bS {input.sam} \
          | samtools view -@ {threads} -b -F 4 - \
          | samtools fastq -@ {threads} \
              -1 {output.out_r1} \
              -2 {output.out_r2} \
              -0 /dev/null -s /dev/null -n
        """

#Assemble mapped reads with SPAdes using k-mer size 99
rule spades_assembly:
    input:
        r1="results/filtered/{sample}_1.fastq.gz",
        r2="results/filtered/{sample}_2.fastq.gz"
    output:
        contigs="results/spades/{sample}/contigs.fasta"
    threads: 8
    shell:
        r"""
        mkdir -p results/spades/{wildcards.sample}

        spades.py \
          -1 {input.r1} \
          -2 {input.r2} \
          -k 99 \
          -t {threads} \
          -o results/spades/{wildcards.sample}
        """

#Calculate number of contigs >1000 bp and total bp
rule assembly_stats:
    input:
        contigs="results/spades/{sample}/contigs.fasta"
    output:
        stats="results/spades/{sample}/assembly_stats.txt"
    shell:
        """
        python scripts/assembly_stats.py \
            -i {input.contigs} \
            -o {output.stats}
        """

#Count raw read pairs (we count reads in R1; paired-end has same count in R2)
rule count_raw_reads:
    input:
        r1="data/raw/{sample}_1.fastq.gz"
    output:
        n="results/counts/{sample}.raw.txt"
    shell:
        """
        mkdir -p results/counts
        python scripts/count_read_pairs.py -i {input.r1} -o {output.n}
        """

#Count filtered read pairs (again count reads in R1)
rule count_filtered_reads:
    input:
        r1="results/filtered/{sample}_1.fastq.gz"
    output:
        n="results/counts/{sample}.filtered.txt"
    shell:
        """
        mkdir -p results/counts
        python scripts/count_read_pairs.py -i {input.r1} -o {output.n}
        """

#Combine raw/filtered read counts + assembly stats into one summary line
rule make_summary_line:
    input:
        raw="results/counts/{sample}.raw.txt",
        filtered="results/counts/{sample}.filtered.txt",
        stats="results/spades/{sample}/assembly_stats.txt",
        strain="results/blast_strain/{sample}.strain.txt"
    output:
        summary="results/summary/{sample}.tsv"
    shell:
        r"""
        mkdir -p results/summary
        python3 -c "raw=open('{input.raw}').read().strip(); filt=open('{input.filtered}').read().strip(); stats=open('{input.stats}').read().strip(); strain=open('{input.strain}').read().strip(); open('{output.summary}','w').write('{wildcards.sample}\t{{}}\t{{}}\t{{}}\t{{}}\n'.format(raw, filt, stats, strain))"
        """

#Convert the per-sample TSV summary into an easy to read report text
rule sample_report:
    input:
        tsv="results/summary/{sample}.tsv"
    output:
        txt="results/report_parts/{sample}.txt"
    shell:
        """
        mkdir -p results/report_parts
        python3 scripts/write_report_one_sample.py -i {input.tsv} -o {output.txt}
        """

#Combine all sample report parts into the final PipelineReport
rule pipeline_report:
    input:
        expand("results/report_parts/{sample}.txt", sample=SAMPLES)
    output:
        "results/PipelineReport.txt"
    shell:
        """
        cat {input} > {output}
        """

#Extract the longest contig from each SPAdes assembly
rule longest_contig:
    input:
        contigs="results/spades/{sample}/contigs.fasta"
    output:
        longest="results/spades/{sample}/longest_contig.fasta"
    shell:
        """
        python3 scripts/longest_contig.py -i {input.contigs} -o {output.longest}
        """

#BLAST the longest contig vs Betaherpesvirinae DB (top 5 hits)
rule blast_top5:
    input:
        query="results/spades/{sample}/longest_contig.fasta",
        db="ref/betaherpes_db/betaherpesvirinae.nsq"
    output:
        tsv="results/blast_strain/{sample}_top5.tsv",
        log="logs/blast_strain/{sample}.log"
    shell:
        """
        mkdir -p results/blast_strain logs/blast_strain
        blastn \
          -query {input.query} \
          -db ref/betaherpes_db/betaherpesvirinae \
          -out {output.tsv} \
          -outfmt "6 qseqid sseqid stitle pident length evalue bitscore" \
          -max_hsps 1 \
          -max_target_seqs 5 \
          2> {output.log}
        """

#Extract the strain name from the top BLAST hit
rule extract_strain:
    input:
        tsv="results/blast_strain/{sample}_top5.tsv"
    output:
        strain="results/blast_strain/{sample}.strain.txt"
    shell:
        """
        python scripts/extract_top_strain.py -i {input.tsv} -o {output.strain}
        """