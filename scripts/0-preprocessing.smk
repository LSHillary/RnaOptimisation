import os
import yaml

# Load samples from the file specified in the configuration
with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

# Define the input files as a list of compressed fastq files
rule all:
    input: "FastQC/Raw/RawMultiQC.html"

# Rule to perform FastQC on raw reads
rule RawFastQC:
    input:
        ForRaw = "0-raw/{sample}_R1_001.fastq.gz",
        RevRaw = "0-raw/{sample}_R2_001.fastq.gz"
    output:
        CheckRawFastQC = "Checks/{sample}_RawFastQC.txt"
    params:
        OutputDirectory = "FastQC/Raw",
        sname = "{sample}"
    threads: 16
    resources:
        mem_mb = 16000,
        partition = "low2",
        cores = 16
    shell:'''
    fastqc {input.ForRaw} -o {params.OutputDirectory} -t {threads} && \
    fastqc {input.RevRaw} -o {params.OutputDirectory} -t {threads} && \
    touch {output.CheckRawFastQC}
    '''

# Rule to run MultiQC on FastQC results
rule RawMultiQC:
    input:
        CheckRawFastQC = expand("Checks/{sample}_RawFastQC.txt", sample = samples),
    output:
        multiqc_report = "FastQC/Raw/RawMultiQC.html"
    params:
        sname = "MultiQC",
        InputDirectory = "FastQC/Raw"
    threads: 2
    resources:
        mem_mb = 8000,
        partition = "low2",
        cores = 2
    shell:'''
        multiqc {params.InputDirectory} --filename {output.multiqc_report}
        '''