import os
import yaml

# Define the rule for filtering contigs and retaining those >= 10Kb in length

# Load samples from the file specified in the configuration
with open(config['samples'], 'r') as f:
    samples = f.read().splitlines()

rule all:
    input: "species_vOTUs.fasta"

rule dRep_vOTUs:
    input:
        input_file="special.fna"
    output:
        vOTU_file="species_vOTUs.fasta"
    threads: 16
    params:
        ani_threshold=95,
        cov_fraction_threshold=85
    shell:
        """
        dRep dereplicate \
            --genomes special.fna \
            --algorithm single_linkage \
            --dereplicateGenomes \
            --genomeInfo special.info.tsv \
            --minANI 95 \
            --minCov 85 \
            --processors 32
        """
