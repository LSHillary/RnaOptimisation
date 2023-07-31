import os
import yaml

# Load the configuration file which will be supplied as argument when running Snakemake
configfile: "pipeline_config.yml"

# Load samples from the file specified in the configuration
with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

# Define the rule for filtering contigs and retaining those >= 10Kb in length
rule all:
    input: expand("Checks/{sample}_genomad_done.txt", sample=samples)

rule filter_10kb:
    input:
        contigsIn = "2-assembly/individual/{sample}.contigs.fa",
        CheckMegahit = "Checks/{sample}_assem_done.txt"
    output:
        contigsOut = "2-assembly/individual/{sample}.contigs.10kb.fa",
        Check10kb = "Checks/{sample}_10kb_done.txt"
    threads: 8
    params:
        sname = "{sample}"
    resources:
        mem_mb = 24000,
        partition = "high2"
    shell:'''
        seqkit seq --min-len 10000 {input.contigsIn} > {output.contigsOut} &&
        touch {output.Check10kb}
        '''

# Define the rule for running GeNomad
rule GeNomad:
    input:
        contigsIn = "2-assembly/individual/{sample}.contigs.10kb.fa",
        Check10kb = "Checks/{sample}_10kb_done.txt"
    output:
        CheckGenomad = "Checks/{sample}_genomad_done.txt"
    threads: 32
    params:
        sname = "{sample}",
        GenomadDB = "/group/jbemersogrp/databases/genomad/genomad_db"
    resources:
        mem_mb = 64000,
        partition = "high2",                                                                                                                                                                                                                                                                                
        time = "24:00:00"
    shell:'''
        genomad end-to-end --cleanup --enable-score-calibration -t {threads} {input.contigsIn} 3-genomad/{params.sname} {params.GenomadDB} &&
        touch {output.CheckGenomad}
        '''
    
