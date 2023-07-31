import os
import yaml

# Load the configuration file which will be supplied as argument when running Snakemake
configfile: "pipeline_config.yml"

# Load samples from the file specified in the configuration
with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

# Define the rule for filtering contigs and retaining those >= 10Kb in length
rule all:
    input: expand("Checks/{sample}_palmscan_done.txt", sample=samples)

# Define the rule for running palmmscan

rule palmscan:
    input:
        contigsIn = "2-assembly/individual/{sample}.contigs.fa",
        CheckMegahit = "Checks/{sample}_assem_done.txt"
    output:
        CheckPalmmscan = "Checks/{sample}_palmscan_done.txt"
    threads: 32
    params:
        sname = "{sample}"
    resources:
        mem_mb = 64000,
        partition = "high2",                                                                                                                                                                                                                                                                                
        time = "24:00:00"
    shell:'''
        source ~/.bashrc
        /group/jbemersogrp/palmscan/bin/palmscan \
            -search_pp {input.contigsIn} \
            -ppout 3-palmscan/{params.sname}_rdrps.faa \
            -ppout_nt 3-palmscan/{params.sname}_rdrps.fna \
            -report 3-palmscan/{params.sname}_pp_report.txt \
            -rdrp -threads 32 && \
        touch {output.CheckPalmmscan}
        '''
    
