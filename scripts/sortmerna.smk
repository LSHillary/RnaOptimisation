# Snakefile for running sortmerna

import os

SAMPLES = [f.replace("_raw_R1.fq.gz", "") for f in os.listdir("0-raw/") if f.endswith("_raw_R1.fq.gz")]

localrules: sortmerna

rule all:
    input:
        expand("Checks/{sample}_sortmerna_done.txt", sample=SAMPLES)

rule sortmerna:
    input:
        ForQC = "1-ProcessedReads/1.1-QC/{sample}_QC_R1.fq.gz",
        RevQC = "1-ProcessedReads/1.1-QC/{sample}_QC_R2.fq.gz",
        CheckQC = "Checks/{sample}_QC_done.txt"
    output:
        CheckSortmeRNA = "Checks/{sample}_sortmerna_done.txt"
    threads: 40
    params:
        mem_mb = 64000,
        sname = "{sample}",
        cores = 40
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
    --out2 \
    --fastx \
    -v --threads 40 && \
    touch {output.CheckSortmeRNA}
    '''
    