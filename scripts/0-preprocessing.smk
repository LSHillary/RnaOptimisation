import os
import yaml

# Load the configuration file which will be supplied as argument when running Snakemake
configfile: "pipeline_config.yml"

# Load samples from the file specified in the configuration
with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

# Define the input files as a list of compressed fastq files
rule all:
    input: "FastQC/Processed/Paired/ProcessedPairedMultiQC.html"

# Rule to merge fastq files from run 1 and run 2
rule MergeRuns:
    input:
        ForRun1 = "0-raw/run_1/{sample}_L004_R1_001.fastq.gz",
        RevRun1 = "0-raw/run_1/{sample}_L004_R2_001.fastq.gz",
        ForRun2 = "0-raw/run_2/{sample}_L004_R1_001.fastq.gz",
        RevRun2 = "0-raw/run_2/{sample}_L004_R2_001.fastq.gz"
    output:
        ForMerge = "0-raw/{sample}_raw_R1.fq.gz",
        RevMerge = "0-raw/{sample}_raw_R2.fq.gz",
        CheckMergeRuns = "Checks/{sample}_MergeRuns.txt"
    params:
        sname = "{sample}"
    threads:1
    resources:
        mem_mb = 24000,
        partition = "high2",
        cores = 1
    message: "Merging runs 1 and 2"
    shell:'''
    cat {input.ForRun1} {input.ForRun2} > {output.ForMerge} && \
    cat {input.RevRun1} {input.RevRun2} > {output.RevMerge} && \
    touch {output.CheckMergeRuns}
    '''

# Rule to perform FastQC on raw reads
rule RawFastQC:
    input:
        ForRaw = "0-raw/{sample}_raw_R1.fq.gz",
        RevRaw = "0-raw/{sample}_raw_R2.fq.gz",
        CheckMergeRuns = "Checks/{sample}_MergeRuns.txt"
    output:
        CheckRawFastQC = "Checks/{sample}_RawFastQC.txt"
    params:
        OutputDirectory = "FastQC/Raw",
        sname = "{sample}"
    threads: 8
    resources:
        mem_mb = 24000,
        partition = "high2",
        cores = 8
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
    threads: 8
    resources:
        mem_mb = 24000,
        partition = "high2",
        cores = 8
    shell:'''
        multiqc {params.InputDirectory} --filename {output.multiqc_report}
        '''
# Adaptor removal and quality filtering on raw reads using bbduk
rule QualityFiltering:
    input:
        ForRaw = "0-raw/{sample}_raw_R1.fq.gz",
        RevRaw = "0-raw/{sample}_raw_R2.fq.gz"
    output:
        ForQC = "1-ProcessedReads/1.1-QC/{sample}_QC_R1.fq.gz",
        RevQC = "1-ProcessedReads/1.1-QC/{sample}_QC_R2.fq.gz",
        SingleQC = "1-ProcessedReads/1.1-QC/{sample}_QC_U.fq.gz",
        CheckQC = "Checks/{sample}_QC_done.txt"
    resources:
        mem_mb = 10000
    threads:8
    params:
        sname = "{sample}"
    message: "QC filtering using bbduk"
    shell:'''
        bbduk.sh in={input.ForRaw} in2={input.RevRaw} \
        ref=adapters,phix ktrim=r k=23 mink=11 hdist=1 tpe tbo \
        qtrim=r trimq=10  maxns=3 maq=3 minlen=50 mlf=0.333 \
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
        CheckQC = "Checks/{sample}_QC_done.txt"
    output:
        CheckProcessedFastQC = "Checks/{sample}_ProcessedFastQC.txt"
    params:
        PairedOutputDirectory = "FastQC/Processed/Paired",
        SingletonsOutputDirectory = "FastQC/Processed/Singletons",
        sname = "{sample}"
    threads: 8
    resources:
        mem_gb = 50
    shell:'''
    fastqc {input.ForProcessed} -o {params.PairedOutputDirectory} -t 8
    fastqc {input.RevProcessed} -o {params.PairedOutputDirectory} -t 8
    fastqc {input.SingleProcessed} -o {params.SingletonsOutputDirectory} -t 8
    touch {output.CheckProcessedFastQC}
    '''
rule ProcessedMultiQC:
    input:
        CheckRawFastQC = expand("Checks/{sample}_ProcessedFastQC.txt", sample = samples)
    output:
        paired_multiqc_report = "FastQC/Processed/Paired/ProcessedPairedMultiQC.html"
    params:
        sname = "MultiQC",
        InputDirectory = "FastQC/Processed/Paired"
    threads: 8
    resources:
        mem_mb = 24000,
        partition = "high2",
        cores = 8
    shell:'''
        multiqc {params.InputDirectory} --filename {output.paired_multiqc_report}
        '''