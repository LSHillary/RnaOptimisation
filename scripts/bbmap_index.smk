import os
import yaml

# Load samples from the file specified in the configuration
with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

rule all:
    input: expand("Checks/bbmap_{sample}_index_done.txt", sample=samples)

rule bbmap_index:
    input:
        contigsIn = "3-genomad/{sample}/{sample}.contigs.10kb_summary/{sample}.contigs.10kb_virus.fna"
    output:
        Index = "4-mapping_bbmap_{sample}_index.idx",
        CheckBBmapIndex = "Checks/bbmap_{sample}_index_done.txt"
    threads: 16
    params:
        sname = "{sample}"
    resources:
        mem_mb = 16000,
        partition = "high2"
    shell:'''
        bbmap.sh ref={input.contigsIn} path=4-mapping && \
        touch {output.CheckBBmapIndex}
        '''