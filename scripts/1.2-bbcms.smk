import os
import yaml

# Load the configuration file which will be supplied as argument when running Snakemake
configfile: "pipeline_config.yml"

# Load samples from the file specified in the configuration
with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

# Define the input files as a list of compressed fastq files
rule all:
    input: expand("Checks/{sample}_bbcms.txt", sample=samples)

# Rule to merge fastq files from run 1 and run 2
rule bbcms:
    input:
        ForQC = "1-ProcessedReads/1.1-QC/{sample}_QC_R1.fq.gz",
        RevQC = "1-ProcessedReads/1.1-QC/{sample}_QC_R2.fq.gz",
        CheckQC = "Checks/{sample}_QC_done.txt"
    output:
        ForBBCMS = "1-ProcessedReads/1.2-bbcms/{sample}_bbcms_R1.fq.gz",
        RevBBCMS = "1-ProcessedReads/1.2-bbcms/{sample}_bbcms_R2.fq.gz",
        ForSingle = "1-ProcessedReads/1.2-bbcms/{sample}_single_R1.fq.gz",
        RevSingle = "1-ProcessedReads/1.2-bbcms/{sample}_single_R2.fq.gz",
        CheckBBCMS = "Checks/{sample}_bbcms.txt"
    params:
        sname = "{sample}"
    threads:12
    resources:
        mem_mb = 24000,
        partition = "high2",
        cores = 12
    message: "Merging runs 1 and 2"
    shell:'''
    bbcms.sh in={input.ForQC} in2={input.RevQC} \
    out={output.ForBBCMS} out2={output.RevBBCMS} \
    outb={output.ForSingle} outb2={output.RevSingle} \
    mincount=2 hcf=0.6 ecc=f k=31 -Xmx24000m && \
    touch {output.CheckBBCMS}
    '''