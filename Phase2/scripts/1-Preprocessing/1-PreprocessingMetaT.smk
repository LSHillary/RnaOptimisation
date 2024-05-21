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
metatranscriptomes = extract_items(data, ['Nucleotides', 'RNA', 'MetaTranscriptomes'])

#### PIPELINE ####

# Rule to check and run all processes
rule all:
    input:
        CheckDnaFastQC = expand("Checks/1.1-RawFastQC_{sample}_MetaT.done", sample=metatranscriptomes),
        CheckDnaMultiQC = "Checks/1.2-MultiQC_MetaT.done",
        CheckPostCoAssemblyCleanup = "Checks/1-PostPreprocessingCleanup.done"

# Rule 1.1 - Run FastQC 
rule RawFastQC:
    input:
        ForRaw = "0-raw/{sample}_L006_R1_001.fastq.gz",
        RevRaw = "0-raw/{sample}_L006_R2_001.fastq.gz",
    output:
        CheckRawFastQC = "Checks/1.1-RawFastQC_{sample}_MetaT.done"
    params:
        OutputDirectory = "1-Preprocessing/1.1-FastQC/Raw/",
        tag = "{sample}_MetaT"
    threads: 8
    resources:
        mem_mb = 32000,
        partition = "high2"
    shell:'''
        mkdir -p {params.OutputDirectory} && \
        fastqc {input.ForRaw} -o {params.OutputDirectory} -t {threads} && \
        fastqc {input.RevRaw} -o {params.OutputDirectory} -t {threads} && \
        touch {output.CheckRawFastQC}
        '''

rule RawMultiQC:
    input:
        CheckRawFastQC = expand("Checks/1.1-RawFastQC_{sample}_MetaT.done", sample = metatranscriptomes)
    output:
        CheckMultiQC = "Checks/1.2-MultiQC_MetaT.done"
    params:
        InputDirectory = "1-Preprocessing/1.1-FastQC/Raw/",
        OutputDirectory = "1-Preprocessing/1.2-MultiQC/Raw/",
        tag = "MetaTRawMultiQC"
    threads: 8
    shell:'''
        mkdir -p {params.OutputDirectory} && \
        multiqc --filename {params.tag} -i {params.tag} -o {params.OutputDirectory} {params.InputDirectory}  && \
        touch {output.CheckMultiQC}
    '''

rule PostPreProcessingCleanup:
    input:
        CheckMultiQC = "Checks/1.2-MultiQC_MetaT.done"
    output:
        CheckPostCoAssemblyCleanup = "Checks/1-PostPreprocessingCleanup.done"
    params:
        tag = "1-PreprocessingCleanup",
        local_folder = "1-Preprocessing/",
        remote_folder = "FARM_archive/RnaOptimisation/Phase2/MetaT/1-Preprocessing/"
    shell:'''
    PASSWORD=$(cat ~/box_password.txt)
    lftp -e "mirror -R --overwrite --only-newer --delete \
    {params.local_folder} {params.remote_folder}; bye" -u lhillary@ucdavis.edu,$PASSWORD ftps://ftp.box.com && \
    rm -r {params.local_folder} && \
    touch {output}
    '''