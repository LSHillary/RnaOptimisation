# QC Snakemake file for CAREERspatial data processing

import os

# Create the SAMPLES list from all the .fastq.gz files in the folder
SAMPLES = [f.replace("_raw_R1.fq.gz", "") for f in os.listdir("0-raw/") if f.endswith("_raw_R1.fq.gz")]

# Define the input files as a list of compressed fastq files
rule all:
    input: "FastQC/Processed/Paired/ProcessedPairedMultiQC.html" # TO EDIT

# Adaptor removal and quality filtering on raw reads using bbduk
rule QualityFiltering:
    input:
        ForProcessed = "1-ProcessedReads/1.1-QC/{sample}_QC_R1.fq.gz",
        RevProcessed = "1-ProcessedReads/1.1-QC/{sample}_QC_R2.fq.gz",
        SingleProcessed = "1-ProcessedReads/1.1-QC/{sample}_QC_U.fq.gz",
        CheckQC = "1-QC/Checks/{sample}_QC_done.txt"
    output:
        ForNorm = "1-ProcessedReads/1.2-Normalisation/{sample}_QC_R1.fq.gz",
        RevNorm = "1-ProcessedReads/1.2-Normalisation/{sample}_QC_R2.fq.gz",
        SingleQC = "1-ProcessedReads/1.2-Normalisation/{sample}_QC_U.fq.gz",
        CheckQC = "1-QC/Checks/{sample}_QC_done.txt"
    resources:
        mem_mb = 10000
    threads:8
    message: "QC filtering using bbduk"
    shell:'''
        bbduk.sh in={input.ForRaw} in2={input.RevRaw} \
        
        out={output.ForQC} out2={output.RevQC} outs={output.SingleQC}
        threads=8 \
        && \
        touch {output.CheckQC}
    '''

rule ProcessedFastQC:
    input:
        ForProcessed = "1-ProcessedReads/1.1-QC/{sample}_QC_R1.fq.gz",
        RevProcessed = "1-ProcessedReads/1.1-QC/{sample}_QC_R2.fq.gz",
        SingleProcessed = "1-ProcessedReads/1.1-QC/{sample}_QC_U.fq.gz",
        CheckQC = "1-QC/Checks/{sample}_QC_done.txt"
    output:
        CheckProcessedFastQC = "Checks/{sample}_ProcessedFastQC.txt"
    params:
        PairedOutputDirectory = "FastQC/Processed/Paired",
        SingletonsOutputDirectory = "FastQC/Processed/Singletons"
    shell:'''
    fastqc {input.ForProcessed} -o {params.PairedOutputDirectory}
    fastqc {input.RevProcessed} -o {params.PairedOutputDirectory}
    fastqc {input.SingleProcessed} -o {params.SingletonsOutputDirectory}
    touch {output.CheckProcessedFastQC}
    '''
rule MultiQC:
    input:
        CheckRawFastQC = expand("Checks/{sample}_ProcessedFastQC.txt", sample = SAMPLES)
    output:
        paired_multiqc_report = "FastQC/Processed/Paired/ProcessedPairedMultiQC.html"

    shell:'''
        multiqc FastQC/Processed/Paired --filename {output.paired_multiqc_report}
        '''
rule BREAK:
    input:
        check = "Checks/RawMultiQC.txt"
    output:
        check = "Checks/done.txt"
    shell:'''
    touch {output.check}'''