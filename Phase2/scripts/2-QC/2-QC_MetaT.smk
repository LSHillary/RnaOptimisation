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
        CheckMultiQC = expand("Checks/2.3-PostFilteringMultiQC.done"),
        #CheckSortmeRNA = expand("Checks/2.4a_sortmerna_{sample}.done", sample=metatranscriptomes),
        CheckRibodetector = expand("Checks/2.4a_ribodetector_{sample}.done", sample = metatranscriptomes),
        CheckPCRDedup = expand("Checks/2.5_Dedup_{sample}.done", sample=metatranscriptomes),
#### 2 - QC Filtering ####

# Rule 2.1 Adaptor removal and quality filtering on raw reads using bbduk
rule QualityFiltering:
    input:
        ForRaw = "0-raw/{sample}_L006_R1_001.fastq.gz",
        RevRaw = "0-raw/{sample}_L006_R2_001.fastq.gz",
    output:
        ForQC="2-QC/2.1-Filtering/{sample}_QC_R1.fq.gz",
        RevQC="2-QC/2.1-Filtering/{sample}_QC_R2.fq.gz",
        CheckQC="Checks/2-QC_{sample}.done"
    resources:
        mem_mb=30000
    threads: 4
    params:
        tag="{sample}"
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
        ForProcessed="2-QC/2.1-Filtering/{sample}_QC_R1.fq.gz",
        RevProcessed="2-QC/2.1-Filtering/{sample}_QC_R2.fq.gz",
        CheckQC="Checks/2-QC_{sample}.done"
    output:
        CheckProcessedFastQC="Checks/FilteredFastQC_{sample}.done"
    params:
        PairedOutputDirectory="2-QC/2.2-FastQC/Filtered/Paired",
        tag="{sample}"
    threads: 4
    resources:
        mem_mb=64000,
        partition="bmh"
    shell:'''
        mkdir -p {params.PairedOutputDirectory} && \
        fastqc {input.ForProcessed} -o {params.PairedOutputDirectory} -t {threads} && \
        fastqc {input.RevProcessed} -o {params.PairedOutputDirectory} -t {threads} && \
        touch {output.CheckProcessedFastQC}
    '''
# Rule 2.3 - MultiQC on filtered reads
rule PostFilteringMultiQC:
    input:
        CheckFilteringFastQC = expand("Checks/FilteredFastQC_{sample}.done", sample=metatranscriptomes)
    output:
        CheckMultiQC = "Checks/2.3-PostFilteringMultiQC.done"
    params:
        InputDirectory = "2-QC/2.2-FastQC/Filtered/Paired",
        OutputDirectory = "2-QC/2.3-MultiQC/Filtered/",
        tag = "FilteredMultiQC"
    threads: 8
    resources:
        mem_mb=128000,
        partition="bmh"
    shell:'''
        mkdir -p {params.OutputDirectory} && \
        multiqc --filename {params.tag} -i {params.tag} -o {params.OutputDirectory} {params.InputDirectory}  && \
        touch {output.CheckMultiQC}
    '''

# Rule 2.4 Error correction using tadpole
rule ErrorCorrection:
    input:
        ForQC="2-QC/2.1-Filtering/{sample}_QC_R1.fq.gz",
        RevQC="2-QC/2.1-Filtering/{sample}_QC_R2.fq.gz"
    output:
        ForEC="2-QC/2.4-ErrorCorrection/{sample}_EC_R1.fq.gz",
        RevEC="2-QC/2.4-ErrorCorrection/{sample}_EC_R2.fq.gz",
        CheckEC="Checks/2.4_EC_{sample}.done"
    resources:
        mem_mb=64000,
        partition = "bmh"
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

rule sortmerna:
    input:
        ForEC="2-QC/2.4-ErrorCorrection/{sample}_EC_R1.fq.gz",
        RevEC="2-QC/2.4-ErrorCorrection/{sample}_EC_R2.fq.gz",
        CheckEC="Checks/2.4_EC_{sample}.done"
    output:
        CheckSortmeRNA = "Checks/2.4a_sortmerna_{sample}.done"
    threads: 32
    resources:
        mem_mb = 131072,
        partition = "bmh",
        time = "2-00:00:00"
    params:
        tag = "{sample}"
    message:
        "Running sortmerna on {input.ForEC} and {input.RevEC}"
    shell:'''
    sortmerna --ref /group/jbemersogrp/databases/sortmerna/smr_v4.3_sensitive_db.fasta \
    --reads {input.ForEC} \
    --reads {input.RevEC} \
    --workdir 2-QC/sortmerna/{params.tag}_wd \
    --idx-dir /home/lhillary/sortmerna/run/idx/ \
    --aligned 2-QC/2.4a-sortmerna/{params.tag}_paired_rrna \
    --other 2-QC/2.4a-sortmerna/{params.tag}_paired_other \
    --paired_in \
    --no-best \
    --num_alignments 1 \
    --out2 \
    --fastx \
    --threads {threads} && \
    touch {output.CheckSortmeRNA}
    '''

rule ribodetector:
    input:
        ForEC="2-QC/2.4-ErrorCorrection/{sample}_EC_R1.fq.gz",
        RevEC="2-QC/2.4-ErrorCorrection/{sample}_EC_R2.fq.gz",
        CheckEC="Checks/2.4_EC_{sample}.done"
    output:
        CheckRibodetector = "Checks/2.4a_ribodetector_{sample}.done",
        ForRD = "2-QC/2.4a-ribodetector/{sample}_RD_R1.fq.gz",
        RevRD = "2-QC/2.4a-ribodetector/{sample}_RD_R2.fq.gz",
        ForRRNA = "2-QC/2.4a-ribodetector/{sample}_RRNA_R1.fq.gz",
        RevRRNA = "2-QC/2.4a-ribodetector/{sample}_RRNA_R2.fq.gz",
        ReadStats = "2-QC/2.4a-ribodetector/{sample}_stats.tsv"
    threads: 32
    resources:
        mem_mb = 32768,
        partition = "bmh",
        time = "2-00:00:00"
    params:
        tag = "{sample}"
    shell:'''
    mkdir -p 2-QC/2.4a-ribodetector && \
    seqkit stats -j16 -T -o {output.ReadStats} {input.ForEC} {input.RevEC} && \
    mean_len=$(awk -F'	' 'NR>1 {{sum+=$7; count++}} END {{if (count > 0) printf "%.0f", sum/count}}' {output.ReadStats}) && \
    micromamba run -n ribodetector ribodetector_cpu \
    -i {input.ForEC} {input.RevEC} \
    -l $mean_len \
    -o {output.ForRD} {output.RevRD} \
    -r {output.ForRRNA} {output.RevRRNA} \
    -e rrna -t 16 --chunk_size 256 && \
    touch {output.CheckRibodetector}
    '''

# Rule 2.5 PCR duplicate removal using clumpify
rule PcrDuplicateRemoval:
    input:
        ForRD = "2-QC/2.4a-ribodetector/{sample}_RD_R1.fq.gz",
        RevRD = "2-QC/2.4a-ribodetector/{sample}_RD_R2.fq.gz",
        CheckRibodetector = "Checks/2.4a_ribodetector_{sample}.done"
    output:
        ForDedup="2-QC/2.5-Deduplication/{sample}_Dedup_R1.fq.gz",
        RevDedup="2-QC/2.5-Deduplication/{sample}_Dedup_R2.fq.gz",
        CheckDedup="Checks/2.5_Dedup_{sample}.done"
    resources:
        mem_mb=64000,
        partition = "bmh"
    threads: 8
    params:
        tag="{sample}",
        OutputDirectory = "2-QC/2.5-Deduplication"
    message: "PCR duplicate removal using clumpify"
    shell:'''
        mkdir -p {params.OutputDirectory} && \
        clumpify.sh \
            in={input.ForRD} \
            in2={input.RevRD} \
            out={output.ForDedup} \
            out2={output.RevDedup} \
            dedupe subs=0 passes=2 \
            threads={threads} \
        && \
        touch {output.CheckDedup}
    '''

rule PostQCCleanup:
    input:
        CheckDedup = expand("Checks/2.5_Dedup_{sample}.done", sample=metatranscriptomes)
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