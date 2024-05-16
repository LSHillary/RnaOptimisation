# Snakemake for running megahit on a set of samples
import os
import yaml

# Load samples from the file specified in the configuration
#with open(config['samples'], 'r') as f:
#    samples = f.read().splitlines()

samples = ['T1_1']

rule all:
    input:
        expand("2-assembly/individual/{sample}.contigs.fa", sample=samples)

rule PE_assembly:
    input:
        CheckQC = "Checks/{sample}_sortmerna_done.txt",
        R1 = "1-ProcessedReads/sortmerna/{sample}_paired_other_fwd.fq.gz",
        R2 = "1-ProcessedReads/sortmerna/{sample}_paired_other_rev.fq.gz"
    output:
        CheckMegahit = "Checks/{sample}_assem_done.txt",
        out_contig = "2-assembly/individual/{sample}.contigs.fa"
    threads: 32
    resources:
        mem_mb = 65536,
        cores = 32,
        partition = "high2",
        time = "8:00:00"
    params:
        tag = "{sample}",
        sname = "{sample}",
        output_folder = "2-assembly/individual",
        output_temp = "2-assembly/megahit_temp"
    message: "paired end assembly on {params.sname}"
    shell:'''
    mkdir -p  2-assembly/megahit_temp/

    # megahit does not allow force overwrite, so each assembly needs to take place in it's own directory. 
    megahit -1 {input.R1} -2 {input.R2} \
    -t 32 --continue --k-min 27 --min-contig-len 300 --continue --presets meta-large \
    --out-dir {params.output_temp}/{wildcards.sample} \
    --out-prefix {wildcards.sample} && \
    mv {params.output_temp}/{wildcards.sample}/{wildcards.sample}.contigs.fa \
    {params.output_folder} && touch {output.CheckMegahit}
    '''
    