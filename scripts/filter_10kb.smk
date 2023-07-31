import os
import yaml

# Define the rule for filtering contigs and retaining those >= 10Kb in length

# Load samples from the file specified in the configuration
with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

rule all:
    input: expand("Checks/{sample}_10kb_done.txt", sample=samples)

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
        module load seqtk
        seqtk seq -L 10000 {input.contigsIn} > {output.contigsOut} && touch {output.Check10kb}
        '''