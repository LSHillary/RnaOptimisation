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
rna_viromes = extract_items(data, ['Nucleotides', 'RNA', 'RnaViromes'])
metatranscriptomes = extract_items(data, ['Nucleotides', 'RNA', 'MetaTranscriptomes'])

DnaRuns = extract_items(data, ['DnaRuns'])
RnaRuns = extract_items(data, ['RnaRuns'])
# Function to map run to runID
def get_DnaRunID(run):
    if run == "run1":
        return "L007"
    elif run == "run2":
        return "L005"
    elif run == "run3":
        return "L001"
    else:
        return "UNKNOWN"  # Handle other cases if necessary

#### PIPELINE ####

# Rule to check and run all processes
rule all:
    input:
        CheckDnaFastQC = expand("Checks/1.1-RawFastQC_{sample}_{run}.done", sample=dna_viromes, run=DnaRuns),
        CheckDnaMultiQC = expand("Checks/1.2-MultiQC_{run}.done", run = DnaRuns),
        CheckPostCoAssemblyCleanup = "Checks/1-PostPreprocessingCleanup.done"

# Rule 1.1 - Run FastQC on merged runs
rule RawFastQC:
    input:
        ForRaw = lambda wildcards: f"0-raw/{wildcards.run}/{wildcards.sample}_{get_DnaRunID(wildcards.run)}_R1_001.fastq.gz",
        RevRaw = lambda wildcards: f"0-raw/{wildcards.run}/{wildcards.sample}_{get_DnaRunID(wildcards.run)}_R2_001.fastq.gz",
    output:
        CheckRawFastQC = "Checks/1.1-RawFastQC_{sample}_{run}.done"
    params:
        OutputDirectory = "1-Preprocessing/1.1-FastQC/Raw/{run}",
        tag = "{sample}"
    threads: 8
    resources:
        mem_mb = 32000,
        partition = "high2"
    message:
        "Running FastQC on {wildcards.sample} for run {wildcards.run}"
    shell:'''
        mkdir -p {params.OutputDirectory} && \
        fastqc {input.ForRaw} -o {params.OutputDirectory} -t {threads} && \
        fastqc {input.RevRaw} -o {params.OutputDirectory} -t {threads} && \
        touch {output.CheckRawFastQC}
        '''

rule RawMultiQC:
    input:
        CheckRawFastQC = expand("Checks/1.1-RawFastQC_{sample}_{run}.done", sample = dna_viromes, run = DnaRuns)
    output:
        CheckMultiQC = "Checks/1.2-MultiQC_{run}.done"
    params:
        InputDirectory = "1-Preprocessing/1.1-FastQC/Raw/{run}",
        OutputDirectory = "1-Preprocessing/1.2-MultiQC/Raw/{run}",
        tag = "{run}RawMultiQC"
    threads: 8
    shell:'''
        mkdir -p {params.OutputDirectory} && \
        multiqc --filename {params.tag} -i {params.tag} -o {params.OutputDirectory} {params.InputDirectory}  && \
        touch {output.CheckMultiQC}
    '''

rule PostPreProcessingCleanup:
    input:
        CheckMultiQC = expand("Checks/1.2-MultiQC_{run}.done", run = DnaRuns)
    output:
        CheckPostCoAssemblyCleanup = "Checks/1-PostPreprocessingCleanup.done"
    params:
        tag = "1-PreprocessingCleanup",
        local_folder = "1-Preprocessing/",
        remote_folder = "FARM_archive/RnaOptimisation/Phase2/Virome/1-Preprocessing/"
    shell:'''
    PASSWORD=$(cat ~/box_password.txt)
    lftp -e "mirror -R --overwrite --only-newer --delete \
    {params.local_folder} {params.remote_folder}; bye" -u lhillary@ucdavis.edu,$PASSWORD ftps://ftp.box.com && \
    rm -r {params.local_folder} && \
    touch {output}
    '''