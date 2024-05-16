#### Section 1 - Preprocessing ####

# Rule 1.1 - Run FastQC on merged runs
rule RawFastQC:
    input:
        ForRaw = "0-raw/run_1/{sample}_L007_R1_001.fastq.gz",
        RevRaw = "0-raw/run_1/{sample}_L007_R2_001.fastq.gz",
    output:
        CheckRawFastQC = "Checks/1.1-RawFastQC_{sample}-run_1.done"
    params:
        OutputDirectory = "1-Preprocessing/FastQC/Raw/run_1",
        tag = "{sample}"
    threads: 8
    resources:
        mem_mb = 32000,
        partition = "high2"
    message:
        "Running FastQC on {wildcards.sample}"
    shell:'''
        mkdir -p {params.OutputDirectory} && \
        fastqc {input.ForRaw} -o {params.OutputDirectory} -t {threads} && \
        fastqc {input.RevRaw} -o {params.OutputDirectory} -t {threads} && \
        touch {output.CheckRawFastQC}
    '''

# # Rule 1.2 - run MultiQC on the FastQC output
# rule RawMultiQC:
#     input:
#         CheckRawFastQC = expand("Checks/1.1-RawFastQC_{sample}.done", sample = samples),
#     output:
#         CheckRawMultiQC = "Checks/1.2-RawMultiQC_All.done",
#         OutputHTML = "RawMultiQC.html"
#     threads: 16
#     resources:
#         mem_mb = 8024
#     params:
#         FastQCDir = "1-Preprocessing/FastQC/Raw/",
#         tag = "Virome"
#     message: "Running MultiQC on Raw Reads"
#     shell:'''
        
#         multiqc {params.FastQCDir}* \
#         --filename {params.tag}_{output.OutputHTML} \
#         --outdir {params.FastQCDir} \
#         --title "Raw MultiQC Report" \
#         && \
#         touch {output.CheckRawMultiQC}
#     '''