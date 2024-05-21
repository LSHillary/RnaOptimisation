import os
import yaml
import glob

# This is the master snakefile that will call all the other snakefiles

# Load the configuration file which will be supplied as argument when running Snakemake
yaml_file_path = "../Config/RnaOptimisationConfigPhase2.yml"

# Open and read the YAML file
with open(yaml_file_path, 'r') as file:
    data = yaml.load(file, Loader=yaml.FullLoader)

def extract_items(yaml_data, path):
    """
    Extracts items from a nested dictionary or list given a list of keys representing the path.
    
    Args:
    yaml_data (dict): The YAML data loaded into a Python dictionary.
    path (list): A list of keys that defines the path to the target items in the nested dictionary.
    
    Returns:
    list: A list of items found at the specified path.
    """
    # Navigate through the nested dictionary using the path
    data_section = yaml_data
    for key in path:
        try:
            data_section = data_section[key]
        except KeyError:
            raise KeyError(f"Key {key} not found in the provided data.")

    # Check if the final section is a dictionary and extract items accordingly
    if isinstance(data_section, dict):
        # Initialize an empty list to collect all items
        all_items = []
        # Loop through each category in the final section and extend the list with its items
        for category in data_section.values():
            all_items.extend(category)
        return all_items
    elif isinstance(data_section, list):
        # Directly return the list if the final node is a list
        return data_section
    else:
        raise ValueError("The data at the specified path is neither a dictionary nor a list.")

# Example usage, adjust paths according to your YAML structure
dna_viromes = extract_items(data, ['Nucleotides', 'DNA', 'DnaViromes'])

runs = extract_items(data, ['Runs'])

# Function to map run to runID
def get_runID(run):
    if run == "run1":
        return "L007"
    elif run == "run2":
        return "L005"
    else:
        return "UNKNOWN"  # Handle other cases if necessary

#### PIPELINE ####

# Rule to check and run all processes
rule all:
    input:
        expand("Checks/4.0-GenomadPrefilter_{sample}.done", sample = dna_viromes),
        CheckCatGenomadResults = "Checks/4.4-CatGenomadResults.done"


rule GeNomadPrefilterindividual:
    input:
        ContigsIn = "3-Assembly/Contigs/{sample}_renamed_contigs.fna",
        CheckRenameContigs = "Checks/3.2-RenameContigs_{sample}.done"
    output:
        TempFasta = "4-virus_identification/genomad/{sample}.fna",
        CheckGenomadPrefilter = "Checks/4.0-GenomadPrefilter_{sample}.done"
    threads: 2
    params:
        tag = "{sample}"
    shell:'''
        seqkit seq -m 1000 {input.ContigsIn} > {output.TempFasta} && \
        touch {output.CheckGenomadPrefilter}
        '''

# Run GeNomad on coassembled contigs
rule GeNomadIndividual:
    input:
        ContigsIn = "4-virus_identification/genomad/{sample}.fna",
        CheckGenomadPrefilter = "Checks/4.0-GenomadPrefilter_{sample}.done"
    output:
        CheckGenomad = "Checks/4.1-Genomad_{sample}.done"
    threads: 12
    params:
        tag = "{sample}",
        GenomadDB = "/group/jbemersogrp/databases/genomad/genomad_db",
        GenomadFolder = "4-virus_identification/genomad/{sample}"
    resources:
        mem_mb = 65536,
        partition = "high2",
        time = "3-00:00:00"
    shell:'''
        micromamba run -n genomad_env genomad end-to-end --cleanup --composition virome --enable-score-calibration -t {threads} {input.ContigsIn} \
        {params.GenomadFolder} {params.GenomadDB} && \
        touch {output.CheckGenomad}
        '''

rule FilterGenomad:
    input:
        CheckGenomad = "Checks/4.1-Genomad_{sample}.done"
    output:
        FilteredGenomadResults = "4-virus_identification/genomad/{sample}/{sample}_FilteredGenomadResults.tsv",
        CheckFilteredGenomadResults = "Checks/4.3-FilterGenomadResults_{sample}.done",
        FilteredContigs = "4-virus_identification/genomad/{sample}_FilteredContigs.fna"
    params:
        tag = "{sample}",
        Contigs = "4-virus_identification/genomad/{sample}/{sample}_summary/{sample}_virus.fna",
        GenomadResults = "4-virus_identification/genomad/{sample}/{sample}_summary/{sample}_virus_summary.tsv"
    shell:'''
        awk -F'\t' 'NR==1 || $2 > 5000 || $11 ~ /Monodnaviria/' {params.GenomadResults} | \
        cut -f1 > {output.FilteredGenomadResults} && \
        seqtk subseq {params.Contigs} {output.FilteredGenomadResults} > {output.FilteredContigs} && \
        touch {output.CheckFilteredGenomadResults}
        '''

rule CatGenomadResults:
    input:
        CheckFilteredGenomadResults = expand("Checks/4.3-FilterGenomadResults_{sample}.done", sample = dna_viromes),
        Contigs = expand("4-virus_identification/genomad/{sample}_FilteredContigs.fna", sample = dna_viromes)
    output:
        CheckCatGenomadResults = "Checks/4.4-CatGenomadResults.done",
        CatGenomadContigs = "4-virus_identification/genomad/individualFilteredViralContigs.fna"
    params:
        tag = "individual"
    shell:'''
        cat {input.Contigs} > {output.CatGenomadContigs} && \
        touch {output.CheckCatGenomadResults}
        '''