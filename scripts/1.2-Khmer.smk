# QC Snakemake file for CAREERspatial data processing

import os

# Create the SAMPLES list from all the .fastq.gz files in the folder
SAMPLES = [f.replace("_raw_R1.fq.gz", "") for f in os.listdir("0-raw/") if f.endswith("_raw_R1.fq.gz")]

# Define the input files as a list of compressed fastq files
rule all:
    input: "FastQC/Khmer/KhmerMultiQC.html" # TO EDIT

rule Interleave:
    input:
        ForProcessed = "1-ProcessedReads/1.1-QC/{sample}_QC_R1.fq.gz",
        RevProcessed = "1-ProcessedReads/1.1-QC/{sample}_QC_R2.fq.gz",
        CheckQC = "Checks/{sample}_QC_done.txt"
    output:
        InterleavedReads = "1-ProcessedReads/1.1-QC/{sample}_I.fq.gz",
        CheckInter = "Checks/{sample}_Inter_done.txt"
    resources:
        mem_mb = 100000
    threads:8
    message: "Interleaving Reads"
    shell:'''
        interleave-reads.py -f -o {output.InterleavedReads} --gzip {input.ForProcessed} {input.RevProcessed} \
        && \
        touch {output.CheckInter}
    '''  
rule khmer:
    input:
        InterleavedReads = "1-ProcessedReads/1.1-QC/{sample}_I.fq.gz",
        CheckInter = "Checks/{sample}_Inter_done.txt"
    output:
        InterKhmer = "1-ProcessedReads/1.2-Khmer/{sample}_Khmer_I.fq.gz",
        CheckKhmer = "Checks/{sample}_Khmer_done.txt"
    resources:
        mem_mb = 600000
    threads:8
    message: "Normalisation using khmer"
    shell:'''
        normalize-by-median.py {input.InterleavedReads} \
        --cutoff=20 --gzip -M 60g \
        --out {output.InterKhmer} \
        && \
        touch {output.CheckKhmer}
    '''  

rule DeInterleave:
    input:
        InterKhmer = "1-ProcessedReads/1.2-Khmer/{sample}_Khmer_I.fq.gz",
        CheckKhmer = "Checks/{sample}_Khmer_done.txt"
    output:
        ForOut = "1-ProcessedReads/1.2-Khmer/{sample}_Khmer_R1.fq.gz",
        RevOut = "1-ProcessedReads/1.2-Khmer/{sample}_Khmer_R2.fq.gz",
        CheckDeInter = "Checks/{sample}_DeInter_done.txt"
    resources:
        mem_mb = 100000
    threads:8
    message: "Interleaving Reads"
    shell:'''
        split-paired-reads.py {input.InterKhmer} -f --gzip -1 {output.ForOut} -2 {output.RevOut} \
        && \
        touch {output.CheckDeInter}
    '''  

rule NormalisedFastQC:
    input:
        ForOut = "1-ProcessedReads/1.2-Khmer/{sample}_Khmer_R1.fq.gz",
        RevOut = "1-ProcessedReads/1.2-Khmer/{sample}_Khmer_R2.fq.gz",
        CheckDeInter = "Checks/{sample}_DeInter_done.txt"
    output:
        CheckNormFastQC = "Checks/{sample}_KhmerFastQC.txt"
    params:
        NormDirectory = "FastQC/Khmer/",
    shell:'''
    fastqc {input.ForOut} -o {params.NormDirectory}
    fastqc {input.RevOut} -o {params.NormDirectory}
    touch {output.CheckNormFastQC}
    '''
rule MultiQC:
    input:
        CheckRawFastQC = expand("Checks/{sample}_KhmerFastQC.txt", sample = SAMPLES)
    output:
        paired_multiqc_report = "FastQC/Khmer/KhmerMultiQC.html"

    shell:'''
        multiqc FastQC/Khmer --filename {output.paired_multiqc_report}
        '''
rule BREAK:
    input:
        check = "Checks/RawMultiQC.txt"
    output:
        check = "Checks/done.txt"
    shell:'''
    touch {output.check}'''
