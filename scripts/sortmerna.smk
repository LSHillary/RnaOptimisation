# Snakefile for running sortmerna

import os
import yaml

# Load samples from the file specified in the configuration
with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

rule all:
    input:
        expand("Checks/{sample}_sortmerna_done.txt", sample=samples)

rule sortmerna:
    input:
        ForQC = "1-ProcessedReads/1.1-QC/{sample}_QC_R1.fq.gz",
        RevQC = "1-ProcessedReads/1.1-QC/{sample}_QC_R2.fq.gz",
        CheckQC = "Checks/{sample}_QC_done.txt"
    output:
        CheckSortmeRNA = "Checks/{sample}_sortmerna_done.txt"
    threads: 32
    resources:
        mem_mb = 131072,
        sname = "{sample}",
        partition = "high2",
        time = "2-00:00:00"
    params:
        sname = "{sample}"
    message:
        "Running sortmerna on {input.ForQC} and {input.RevQC}"
    shell:'''
    sortmerna --ref /group/jbemersogrp/databases/sortmerna/smr_v4.3_sensitive_db.fasta \
    --reads {input.ForQC} \
    --reads {input.RevQC} \
    --workdir 1-ProcessedReads/sortmerna/{params.sname}_wd \
    --idx-dir /home/lhillary/sortmerna/run/idx/ \
    --aligned 1-ProcessedReads/sortmerna/{params.sname}_paired_rrna \
    --other 1-ProcessedReads/sortmerna/{params.sname}_paired_other \
    --paired_in \
    --num_alignments 1 \
    --out2 \
    --fastx \
    -v --threads 40 && \
    touch {output.CheckSortmeRNA}
    '''
    