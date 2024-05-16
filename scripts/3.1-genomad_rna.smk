import os
import yaml

# Load samples from the file specified in the configuration
with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

# Define the rule for filtering contigs and retaining those >= 10Kb in length
rule all:
    input:
        GeNomad_Results = expand("Checks/{sample}_genomad_done.txt", sample=samples),
        PalmScan_Results = expand("Checks/{sample}_palmscan_done.txt", sample=samples)

# Define the rule for running GeNomad
rule GeNomad:
    input:
        contigsIn = "2-assembly/individual/{sample}.contigs.fa",
        Check = "Checks/{sample}_assem_done.txt"
    output:
        CheckGenomad = "Checks/{sample}_genomad_done.txt"
    threads: 32
    params:
        tag = "{sample}",
        sname = "{sample}",
        GenomadDB = "/group/jbemersogrp/databases/genomad/genomad_db"
    resources:
        mem_mb = 64000,
        partition = "high2",                                                                                                                                                                                                                                                                                
        time = "24:00:00"
    shell:'''
        mkdir -p 3-genomad/{params.sname} && \
        genomad end-to-end --cleanup --enable-score-calibration -t {threads} {input.contigsIn} 3-genomad/{params.sname} {params.GenomadDB} &&
        touch {output.CheckGenomad}
        '''


# Define the rule for running palmmscan

rule palmscan:
    input:
        contigsIn = "2-assembly/individual/{sample}.contigs.fa",
        CheckMegahit = "Checks/{sample}_assem_done.txt"
    output:
        CheckPalmmscan = "Checks/{sample}_palmscan_done.txt"
    threads: 32
    params:
        sname = "{sample}",
        tag = "{sample}"
    resources:
        mem_mb = 64000,
        partition = "high2",                                                                                                                                                                                                                                                                                
        time = "24:00:00"
    shell:'''
        source ~/.bashrc
        mkdir -p 3-palmscan && \
        /group/jbemersogrp/palmscan/bin/palmscan \
            -search_pp {input.contigsIn} \
            -ppout 3-palmscan/{params.sname}_rdrps.faa \
            -ppout_nt 3-palmscan/{params.sname}_rdrps.fna \
            -report 3-palmscan/{params.sname}_pp_report.txt \
            -rdrp -threads 32 && \
        touch {output.CheckPalmmscan}
        '''
    
