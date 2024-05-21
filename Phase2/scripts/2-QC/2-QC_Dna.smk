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
        CheckMultiQC = expand("Checks/2.3-PostFilteringMultiQC_{run}.done", run = DnaRuns),
        CheckPCRDedup = expand("Checks/2.5_Dedup_{sample}_{run}.done", sample = dna_viromes, run = DnaRuns),
#### 2 - QC Filtering ####

# Rule 2.1 Adaptor removal and quality filtering on raw reads using bbduk
rule QualityFiltering:
    input:
        ForRaw = lambda wildcards: f"0-raw/{wildcards.run}/{wildcards.sample}_{get_DnaRunID(wildcards.run)}_R1_001.fastq.gz",
        RevRaw = lambda wildcards: f"0-raw/{wildcards.run}/{wildcards.sample}_{get_DnaRunID(wildcards.run)}_R2_001.fastq.gz",
    output:
        ForQC="2-QC/2.1-Filtering/{sample}_{run}_QC_R1.fq.gz",
        RevQC="2-QC/2.1-Filtering/{sample}_{run}_QC_R2.fq.gz",
        CheckQC="Checks/2-QC_{sample}_{run}.done"
    resources:
        mem_mb=10000
    threads: 4
    params:
        tag="{sample}_{run}"
    message: "QC filtering using bbduk"
    shell:'''
        mkdir -p 2-QC/2.1-Filtering && \
        bbduk.sh in={input.ForRaw} in2={input.RevRaw} \
        ref=adapters,phix ktrim=r k=23 mink=10 rcomp=t hdist=1 tpe tbo threads={threads} \
        qtrim=r trimq=10 maxns=3 maq=3 minlen=50 mlf=0.333 mingc=0.05 maxgc=0.95 \
        out={output.ForQC} out2={output.RevQC} \
        && \
        touch {output.CheckQC}
    '''

# Rule 2.2 Run FastQC on processed reads
rule FilteringFastQC:
    input:
        ForProcessed="2-QC/2.1-Filtering/{sample}_{run}_QC_R1.fq.gz",
        RevProcessed="2-QC/2.1-Filtering/{sample}_{run}_QC_R2.fq.gz",
        CheckQC="Checks/2-QC_{sample}_{run}.done"
    output:
        CheckProcessedFastQC="Checks/FilteredFastQC_{sample}_{run}.done"
    params:
        PairedOutputDirectory="2-QC/2.2-FastQC/Filtered/Paired",
        tag="{sample}"
    threads: 4
    resources:
        mem_gb=30,
        partition="high2"
    shell:'''
        mkdir -p {params.PairedOutputDirectory} && \
        fastqc {input.ForProcessed} -o {params.PairedOutputDirectory} -t {threads} && \
        fastqc {input.RevProcessed} -o {params.PairedOutputDirectory} -t {threads} && \
        touch {output.CheckProcessedFastQC}
    '''
# Rule 2.3 - MultiQC on filtered reads
rule PostFilteringMultiQC:
    input:
        CheckFilteringFastQC = expand("Checks/FilteredFastQC_{sample}_{run}.done", sample = dna_viromes, run = DnaRuns)
    output:
        CheckMultiQC = "Checks/2.3-PostFilteringMultiQC_{run}.done"
    params:
        InputDirectory = "2-QC/2.2-FastQC/Filtered/Paired",
        OutputDirectory = "2-QC/2.3-MultiQC/Filtered/{run}",
        tag = "{run}_FilteredMultiQC"
    threads: 8
    shell:'''
        mkdir -p {params.OutputDirectory} && \
        multiqc --filename {params.tag} -i {params.tag} -o {params.OutputDirectory} {params.InputDirectory}  && \
        rm -r 2-QC/2.2-MultiQC/Filtered && \
        touch {output.CheckMultiQC}
    '''

# Rule 2.4 Error correction using tadpole
rule ErrorCorrection:
    input:
        ForQC="2-QC/2.1-Filtering/{sample}_{run}_QC_R1.fq.gz",
        RevQC="2-QC/2.1-Filtering/{sample}_{run}_QC_R2.fq.gz"
    output:
        ForEC="2-QC/2.4-ErrorCorrection/{sample}_{run}_EC_R1.fq.gz",
        RevEC="2-QC/2.4-ErrorCorrection/{sample}_{run}_EC_R2.fq.gz",
        CheckEC="Checks/2.4_EC_{sample}_{run}.done"
    resources:
        mem_mb=32000,
        partition = "high2"
    threads: 8
    params:
        tag="{sample}",
        OutputDirectory = "2-QC/2.4-ErrorCorrection"
    message: "Error correction using tadpole"
    shell:'''
    mkdir -p {params.OutputDirectory} && \
    tadpole.sh \
        in={input.ForQC} \
        in2={input.RevQC} \
        out={output.ForEC} \
        out2={output.RevEC} \
        mode=correct ecc=t prefilter=1 \
        threads={threads} \
    && \
    touch {output.CheckEC}
    '''

# Rule 2.5 PCR duplicate removal using clumpify
rule PcrDuplicateRemoval:
    input:
        ForEC="2-QC/2.4-ErrorCorrection/{sample}_{run}_EC_R1.fq.gz",
        RevEC="2-QC/2.4-ErrorCorrection/{sample}_{run}_EC_R2.fq.gz",
        CheckEC="Checks/2.4_EC_{sample}_{run}.done"
    output:
        ForDedup="2-QC/2.5-Deduplication/{sample}_{run}_Dedup_R1.fq.gz",
        RevDedup="2-QC/2.5-Deduplication/{sample}_{run}_Dedup_R2.fq.gz",
        CheckDedup="Checks/2.5_Dedup_{sample}_{run}.done"
    resources:
        mem_mb=24000,
        partition = "high2"
    threads: 8
    params:
        tag="{sample}",
        OutputDirectory = "2-QC/2.5-Deduplication"
    message: "PCR duplicate removal using clumpify"
    shell:'''
        mkdir -p {params.OutputDirectory} && \
        clumpify.sh \
            in={input.ForEC} \
            in2={input.RevEC} \
            out={output.ForDedup} \
            out2={output.RevDedup} \
            dedupe subs=0 passes=2 \
            threads={threads} \
        && \
        touch {output.CheckDedup}
    '''

rule PostQCCleanup:
    input:
        CheckDedup = expand("Checks/2.5_Dedup_{sample}_{run}.done", sample = dna_viromes, run = DnaRuns)
    output:
        CheckQCCleanup = "Checks/2-PostQCCleanup.done"
    params:
        tag = "2-QCCleanup",
        local_folder = "2-QC/2.2-MultiQC/",
        remote_folder = "FARM_archive/TestData/2-QC/MultiQC"
    shell:'''
    PASSWORD=$(cat ~/box_password.txt)
    lftp -e "mirror -R --overwrite --only-newer --delete \
    {params.local_folder} {params.remote_folder}; bye" -u lhillary@ucdavis.edu,$PASSWORD ftps://ftp.box.com && \
    rm -r {params.local_folder} && \
    touch {output.CheckQCCleanup}
    '''