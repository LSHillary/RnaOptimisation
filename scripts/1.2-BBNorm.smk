# QC Snakemake file for CAREERspatial data processing

import os

# Create the SAMPLES list from all the .fastq.gz files in the folder
SAMPLES = [f.replace("_raw_R1.fq.gz", "") for f in os.listdir("0-raw/") if f.endswith("_raw_R1.fq.gz")]

# Define the input files as a list of compressed fastq files
rule all:
    input: "FastQC/BBNorm/BBNormMultiQC.html" # TO EDIT

# Adaptor removal and quality filtering on raw reads using bbduk
rule BBNorm:
    input:
        ForProcessed = "1-ProcessedReads/1.1-QC/{sample}_QC_R1.fq.gz",
        RevProcessed = "1-ProcessedReads/1.1-QC/{sample}_QC_R2.fq.gz",
        CheckQC = "Checks/{sample}_QC_done.txt"
    output:
        ForNorm = "1-ProcessedReads/1.2-BBNorm/{sample}_BBNorm_R1.fq.gz",
        RevNorm = "1-ProcessedReads/1.2-BBNorm/{sample}_BBNorm_R2.fq.gz",
        CheckNorm = "Checks/{sample}_BBNorm_done.txt"
    resources:
        mem_mb = 60000,
        cores = 40,
        partition = "high2"
    threads:40
    params:
        sname = "{sample}"
    message: "Normalisation using bbnorm"
    shell:'''
        bbnorm.sh in={input.ForProcessed} in2={input.RevProcessed} \
        target=100 min=1 \
        out={output.ForNorm} out2={output.RevNorm} \
        threads=40 -Xmx6000 \
        && \
        touch {output.CheckNorm}
    '''

rule BBNormFastQC:
    input:
        ForNorm = "1-ProcessedReads/1.2-BBNorm/{sample}_BBNorm_R1.fq.gz",
        RevNorm = "1-ProcessedReads/1.2-BBNorm/{sample}_BBNorm_R2.fq.gz",
        CheckQC = "Checks/{sample}_BBNorm_done.txt"
    output:
        CheckNormFastQC = "Checks/{sample}_BBNormFastQC.txt"
    params:
        NormDirectory = "FastQC/BBNorm/",
        sname = "{sample}"
    shell:'''
    fastqc {input.ForNorm} -o {params.NormDirectory}
    fastqc {input.RevNorm} -o {params.NormDirectory}
    touch {output.CheckNormFastQC}
    '''
rule MultiQC:
    input:
        CheckRawFastQC = expand("Checks/{sample}_BBNormFastQC.txt", sample = SAMPLES)
    output:
        BBNorm_multiqc_report = "FastQC/BBNorm/BBNormMultiQC.html"
    params:
        sname = "MultiQC"
    shell:'''
        multiqc FastQC/BBNorm --filename {output.BBNorm_multiqc_report}
        '''