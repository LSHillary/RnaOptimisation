import os
import yaml

# Load samples from the file specified in the configuration
with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

rule all:
    input: expand("Checks/{sample}_genomad_done.txt", sample=samples)

# Define the rule for running GeNomad
rule GeNomad:
    input:
        contigsIn = "2-assembly/individual/{sample}.contigs.10kb.fa",
        Check10kb = "Checks/{sample}_10kb_done.txt"
    output:
        CheckGenomad = "Checks/{sample}_genomad_done.txt"
    threads: 16
    params:
        sname = "{sample}",
        GenomadDB = "/group/jbemersogrp/databases/genomad/genomad_db"
    resources:
        mem_mb = 32000,
        partition = "low2"
    shell:'''
        mkdir -p 3-genomad/{params.sname} &&
        genomad end-to-end --cleanup --enable-score-calibration -t {threads} {input.contigsIn} 3-genomad/{params.sname} {params.GenomadDB} &&
        touch {output.CheckGenomad}
        '''
    
