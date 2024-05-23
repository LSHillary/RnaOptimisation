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
        #CheckAssembly = expand("Checks/3-PostQCCleanup_{sample}.done", sample = metatranscriptomes),
        CheckRenameAssemblyContigs = expand("Checks/3.2-RenameContigs_{sample}.done", sample = metatranscriptomes),
        CheckGenomad = expand("Checks/4.1-Genomad_{sample}.done", sample = metatranscriptomes),

rule ErrorCorrection:
    input:
        ForQC="2-QC/2.1-Filtering/{sample}_QC_R1.fq.gz",
        RevQC="2-QC/2.1-Filtering/{sample}_QC_R2.fq.gz"
    output:
        ForEC="2-QC/2.4-ErrorCorrection/{sample}_EC_R1.fq.gz",
        RevEC="2-QC/2.4-ErrorCorrection/{sample}_EC_R2.fq.gz",
        CheckEC="Checks/2.4_EC_{sample}.done"
    resources:
        mem_mb=100000,
        partition = "bmh"
    threads: 16
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

# Assembly of all samples by bucket using megahit
rule Assembly_megahit:
    input:
        For1="2-QC/2.5-Deduplication/{sample}_Dedup_R1.fq.gz",
        Rev1="2-QC/2.5-Deduplication/{sample}_Dedup_R2.fq.gz",
    output:
        CheckMegahitAssembly="Checks/3.1-Assembly_{sample}.done"
    params:
        tag="{sample}",
        output_folder="3-Assembly/Contigs",
        output_temp="3-Assembly/megahit_temp",
    threads: 16
    resources:
        mem_mb=48000,
        partition="bmh",
        time="2-00:00:00"
    shell:'''
        mkdir -p {params.output_temp} && \
        megahit -1 {input.For1} -2 {input.Rev1} \
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

rule GeNomadIndividual:
    input:
        ContigsIn = "3-Assembly/Contigs/{sample}_renamed_contigs.fna",
        CheckRenameAssemblyContigs = "Checks/3.2-RenameContigs_{sample}.done"
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
        CheckRenameAssemblyContigs = expand("Checks/3.4-RenameContigs_{sample}.done", sample=metatranscriptomes),
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