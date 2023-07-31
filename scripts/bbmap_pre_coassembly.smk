import os
import yaml

# Load samples from the file specified in the configuration
with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

rule all:
    input: expand("Checks/{sample}_bbmap_combine_done.txt", sample=samples)

# Define the rule for running minimap2
rule bbmap_align:
    input:
        ForQC = "1-ProcessedReads/1.1-QC/{sample}_QC_R1.fq.gz",
        RevQC = "1-ProcessedReads/1.1-QC/{sample}_QC_R2.fq.gz"
    output:
        CheckBBmapAlign = "Checks/{sample}_bbmap_combine_done.txt",
        ForUnmapped = "2-assembly/coassembly/{sample}unmapped_R1.fq.gz",
        RevUnmapped = "2-assembly/coassembly/{sample}unmapped_R2.fq.gz"
    threads: 16
    params:
        sname = "{sample}",
        index_path = "2-assembly/individual/PreCoassembly"
    resources:
        mem_mb = 32000,
        partition = "high2"
    shell:'''
        bbmap.sh in1={input.ForQC} in2={input.RevQC} \
        minid=0.9 vslow=t path={params.index_path} \
        outu={output.ForUnmapped} outu2={output.RevUnmapped} && \
        touch {output.CheckBBmapAlign}
        '''