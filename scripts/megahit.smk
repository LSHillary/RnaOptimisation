# Snakemake for running megahit on a set of samples

import os

# Load the configuration file which will be supplied as argument when running Snakemake
configfile: "pipeline_config.yml"

# Load samples from the file specified in the configuration
with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

rule all:
    input:
        expand("2-assembly/individual/{sample}.contigs.fa", sample=samples)

rule PE_assembly:
    input:
        CheckQC = "Checks/{sample}_QC_done.txt",
        R1 = "1-ProcessedReads/1.1-QC/{sample}_QC_R1.fq.gz",
        R2 = "1-ProcessedReads/1.1-QC/{sample}_QC_R2.fq.gz"
    output:
        CheckMegahit = "Checks/{sample}_assem_done.txt",
        out_contig = "2-assembly/individual/{sample}.contigs.fa"
    threads: 16
    resources:
        mem_mb = 32000,
        cores = 16,
        partition = bmh,
        time = 5-24:00:00
    params:
        sname = "{sample}",
        output_folder = "2-assembly/individual",
        output_temp = "2-assembly/megahit_temp"
    message: "paired end assembly on {params.sname}"
    shell:'''
    mkdir -p  2-assembly/megahit_temp/

    # megahit does not allow force overwrite, so each assembly needs to take place in it's own directory. 
    megahit -1 {input.R1} -2 {input.R2} \
    -t 16 --continue --k-min 27 --min-contig-len 1000 --presets meta-large \
    --out-dir {params.output_temp}/{wildcards.sample} \
    --out-prefix {wildcards.sample} && \
    mv {params.output_temp}/{wildcards.sample}/{wildcards.sample}.contigs.fa \
    {params.output_folder} && touch {output.CheckMegahit}
    '''
    