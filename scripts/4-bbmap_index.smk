import os
import yaml

# Load samples from the file specified in the configuration
with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

rule all:
    input: expand("Checks/bbmap_{sample}_mapping_done.txt", sample=samples)

rule bbmap_index:
    input:
        contigsIn = "3-genomad/{sample}/{sample}.contigs_summary/{sample}.contigs_virus.fna"
    output:
        OutputFolder = "4-mapping/{sample}/",
        CheckBBmapIndex = "Checks/bbmap_{sample}_index_done.txt"
    threads: 16
    params:
        sname = "{sample}"
    resources:
        mem_mb = 16000,
        partition = "high2"
    shell:'''
        bbmap.sh ref={input.contigsIn} path=4-mapping/{wildcards.sample}/ && \
        touch {output.CheckBBmapIndex}
        '''

rule bbmap_map:
    input:
        contigsIn = "3-genomad/{sample}/{sample}.contigs_summary/{sample}.contigs_virus.fna",
        FwdReads = "1-ProcessedReads/sortmerna/{sample}_paired_other_fwd.fq.gz",
        RevReads = "1-ProcessedReads/sortmerna/{sample}_paired_other_rev.fq.gz"
    output:
        sam = "4-mapping/{sample}/{sample}_mapped.sam",
        covstats = "4-mapping/{sample}/{sample}_coverage_stats.txt",
        rpkm = "4-mapping/{sample}/{sample}_rpkm.txt",
        CheckBBmapMapping = "Checks/bbmap_{sample}_mapping_done.txt"
    threads: 16
    params:
        sname = "{sample}"
    resources:
        mem_mb = 16000,
        partition = "high2"
    shell:'''
        bbmap.sh in1={input.FwdReads} in2={input.RevReads} ref={input.contigsIn} \
        out={output.sam} covstats={output.covstats} rpkm={output.rpkm} minid=0.9 threads={threads} && \
        touch {output.CheckBBmapMapping}
        '''
