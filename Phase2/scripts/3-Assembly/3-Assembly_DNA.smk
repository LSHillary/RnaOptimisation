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

#### PIPELINE ####

# Rule to check and run all processes
rule all:
    input:
        #CheckAssembly = expand("Checks/3-PostQCCleanup_{sample}.done", sample = dna_viromes),
        CheckRenameAssemblyContigs = expand("Checks/3.2-RenameContigs_{sample}.done", sample = dna_viromes),

# Assembly of all samples by bucket using megahit
rule Assembly_megahit:
    input:
        For1="2-QC/2.5-Deduplication/{sample}_run1_Dedup_R1.fq.gz",
        Rev1="2-QC/2.5-Deduplication/{sample}_run1_Dedup_R2.fq.gz",
        For2="2-QC/2.5-Deduplication/{sample}_run2_Dedup_R1.fq.gz",
        Rev2="2-QC/2.5-Deduplication/{sample}_run2_Dedup_R2.fq.gz",
        For3="2-QC/2.5-Deduplication/{sample}_run3_Dedup_R1.fq.gz",
        Rev3="2-QC/2.5-Deduplication/{sample}_run3_Dedup_R2.fq.gz",
    output:
        CheckMegahitAssembly="Checks/3.1-Assembly_{sample}.done"
    params:
        tag="{sample}",
        output_folder="3-Assembly/Contigs",
        output_temp="3-Assembly/megahit_temp",
    threads: 16
    resources:
        mem_mb=48000,
        partition="high2",
        time="2-00:00:00"
    shell:'''
        mkdir -p {params.output_temp} && \
        megahit -1 {input.For1},{input.For2},{input.For3} -2 {input.Rev1},{input.Rev2},{input.Rev3} \
        -t {threads} --continue --k-min 27 --presets meta-large \
        --out-dir {params.output_temp}/{params.tag} \
        --out-prefix {params.tag} && touch {output.CheckMegahitAssembly}
        '''


# Rename coassembled contigs
rule RenameAssemblyContigs:
    input:
        CheckMegahitAssembly = "Checks/3.1-Assembly_{sample}.done"
    output:
        renamed_contigs = "3-Assembly/Contigs/{sample}_renamed_contigs.fna",
        CheckRenameAssemblyContigs = "Checks/3.2-RenameContigs_{sample}.done"
    params:
        tag = "{sample}",
        contigs = "3-Assembly/megahit_temp/{sample}/{sample}.contigs.fa"
    shell:'''
    awk '/^>/{{print ">" "{params.tag}" "_contig_" ++i; next}}{{print}}' < {params.contigs} > {output.renamed_contigs} && \
    touch {output.CheckRenameAssemblyContigs}
    '''
# Clean up and archive all intermediate files not needed in future steps
rule PostAssemblyCleanup:
    input:
        CheckRenameAssemblyContigs = "Checks/3.2-RenameContigs_{sample}.done"
    output:
        CheckAssemblyCleanup = "Checks/3-PostQCCleanup_{sample}.done"
    params:
        tag = "{sample}",
        local_Assembly = "3-Assembly/megahit_temp/{sample}",
        remote_Assembly = "FARM_archive/TestData/3-Assembly/megahit_temp/{sample}",
    shell:'''
    PASSWORD=$(cat ~/box_password.txt)
    lftp -e "mirror -R --overwrite --only-newer --delete \
    {params.local_Assembly} {params.remote_Assembly}; bye" -u lhillary@ucdavis.edu,$PASSWORD ftps://ftp.box.com && \
    rm -r {params.local_Assembly} && \
    touch {output.CheckAssemblyCleanup}
    '''

# Rule to push assembly reads to Box
rule PushAssemblyReadsToBox:
    input:
        CheckRenameAssemblyContigs = expand("Checks/3.4-RenameContigs_{sample}.done", sample=dna_viromes),
    output:
        CheckAssemblyReadsCleanup = "Checks/3-ReadsUpload.done"
    params:
        tag = "PushReadsToBox",
        local_reads_dir = "2-QC/2.5-Deduplication/",
        local_filtered_reads_dir = "2-QC/2.1-Filtering",
        remote_reads_dir = "FARM_archive/TestData/2-QC/2.5-Deduplication/",
        remote_filtered_reads_dir = "FARM_archive/TestData/2-QC/2.1-Filtering/"
    shell:
        '''
        PASSWORD=$(cat ~/box_password.txt)
        #lftp -e "mirror -R --overwrite --only-newer --delete \
        #{params.local_reads_dir} {params.remote_reads_dir}; bye" -u lhillary@ucdavis.edu,$PASSWORD ftps://ftp.box.com && \
        #rm -r {params.local_reads_dir} && \
        && \
        lftp -e "mirror -R --overwrite --only-newer --delete \
        {params.local_filtered_reads_dir} {params.remote_filtered_reads_dir}; bye" -u lhillary@ucdavis.edu,$PASSWORD ftps://ftp.box.com && \
        #rm -r {params.local_reads_dir} && \
        rm -r {params.local_filtered_reads_dir} && \
        touch {output.CheckAssemblyReadsCleanup}
        '''