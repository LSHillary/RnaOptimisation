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
rna_viromes = extract_items(data, ['Nucleotides', 'RNA', 'RnaViromes'])

#### PIPELINE ####

# Rule to check and run all processes
rule all:
    input:
        #CheckAssembly = expand("Checks/3-PostQCCleanup_{sample}.done", sample = dna_viromes),
        #CheckRenameAssemblyContigs = expand("Checks/3.2-RenameContigs_{sample}.done", sample = rna_viromes),
        CheckGenomad = expand("Checks/4.1-Genomad_{sample}.done", sample = rna_viromes),
        CheckGenomadRerun = expand("Checks/4.4_GenomadRerun_{sample}.done", sample = rna_viromes),

# Assembly of all samples by bucket using megahit
rule Assembly_megahit:
    input:
        For1="2-QC/2.5-Deduplication/{sample}_run1_Dedup_R1.fq.gz",
        Rev1="2-QC/2.5-Deduplication/{sample}_run1_Dedup_R2.fq.gz",
        #For2="2-QC/2.5-Deduplication/{sample}_run2_Dedup_R1.fq.gz",
        #Rev2="2-QC/2.5-Deduplication/{sample}_run2_Dedup_R2.fq.gz",
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
        mkdir -p {params.GenomadFolder} && \
        micromamba run -n genomad_env genomad end-to-end --cleanup --enable-score-calibration -t {threads} {input.ContigsIn} \
        {params.GenomadFolder} {params.GenomadDB} && \
        touch {output.CheckGenomad}
        '''

#### 4.2 - Contig Extension ####
# 4.2.1 - Mapping
rule CE_Mapping:
    input:
        FwdReadsRun1="2-QC/2.4a-ribodetector/{sample}_run1_RD_R1.fq.gz",
        RevReadsRun1="2-QC/2.4a-ribodetector/{sample}_run1_RD_R2.fq.gz",
        #FwdReadsRun2="2-QC/2.4a-ribodetector/{sample}_run2_EC_R1.fq.gz",
        #RevReadsRun2="2-QC/2.4a-ribodetector/{sample}_run2_EC_R2.fq.gz",
        Contigs = "3-Assembly/Contigs/{sample}_renamed_contigs.fna",
        CheckGenomad = "Checks/4.1-Genomad_{sample}.done"
    output:
        Bam=temp("4-virus_identification/4.3-ContigExtension/{sample}_mapped.bam"),
        FwdReadsAll=temp("4-virus_identification/4.3-ContigExtension/{sample}_mapped_R1.fq.gz"),
        RevReadsAll=temp("4-virus_identification/4.3-ContigExtension/{sample}_mapped_R2.fq.gz"),
        CheckMinimap2Mapping="Checks/4.3-ContigExtension_{sample}.done"
    threads: 16
    params:
        tag = "{sample}"
    resources:
        mem_mb = 128000,
        partition = "high2",
        time = "3-00:00:00"
    shell:'''
        cat {input.FwdReadsRun1}  > {output.FwdReadsAll} && \
        cat {input.RevReadsRun1}  > {output.RevReadsAll} && \
        minimap2 -axsr -t {threads} {input.Contigs} {output.FwdReadsAll} {output.RevReadsAll} | samtools view -u | samtools sort -o {output.Bam} \
        && \
        touch {output.CheckMinimap2Mapping}
        '''
rule CE_CoverM:
    input:
        Bam="4-virus_identification/4.3-ContigExtension/{sample}_mapped.bam",
        CheckMinimap2Mapping="Checks/4.3-ContigExtension_{sample}.done"
    output:
        Coverage = "4-virus_identification/4.3-ContigExtension/{sample}_coverage.txt",
        CheckCoverM = "Checks/4.3-ContigExtension_{sample}_CoverM.done"
    threads: 16
    params:
        tag = "{sample}",
        CoverMPath = "4-virus_identification/4.3-ContigExtension/",
        min_coverage = 75
    resources:
        mem_mb = 16000,
        partition = "high2",
        time = "19:00:00"
    shell:'''
        coverm contig -t {threads} --methods metabat --min-read-percent-identity 90 \
        --bam-files {input.Bam} --output-format dense -o {output.Coverage} && \
        touch {output.CheckCoverM}
        '''
rule CE_filtering:
    input:
        Fasta = "3-Assembly/Contigs/{sample}_renamed_contigs.fna",
        Coverage = "4-virus_identification/4.3-ContigExtension/{sample}_coverage.txt",
        CheckCoverM = "Checks/4.3-ContigExtension_{sample}_CoverM.done"
    output:
        Query = temp("4-virus_identification/4.3-ContigExtension/{sample}_query.fa"),
        QueryContigNames = temp("4-virus_identification/4.3-ContigExtension/{sample}_query_contig_names.txt"),
        CheckCE_filtering = "Checks/4.3-ContigExtension_{sample}_filtering.done",
        CoverageDoc = "4-virus_identification/4.3-ContigExtension/{sample}_coverage_new.txt"
    threads: 4
    resources:
        mem_mb = 8000,
        partition = "high2"
    params:
        tag = "{sample}",
        GenomadFasta = "4-virus_identification/genomad/{sample}/{sample}_renamed_contigs_summary/{sample}_renamed_contigs_virus.fna"
    shell:'''
        grep '^>' {params.GenomadFasta} | cut -d'|' -f1 | sed 's/>//' > {output.QueryContigNames} && \
        seqtk subseq {input.Fasta} {output.QueryContigNames} > {output.Query} && 
        python ../scripts/coverage.transfer.py -i {input.Coverage} -o {output.CoverageDoc} && \
        touch {output.CheckCE_filtering}
        '''
rule ContigExtension:
    input:
        Fasta = "3-Assembly/Contigs/{sample}_renamed_contigs.fna",
        Query = "4-virus_identification/4.3-ContigExtension/{sample}_query.fa",
        CoverageDoc = "4-virus_identification/4.3-ContigExtension/{sample}_coverage_new.txt",
        Bam="4-virus_identification/4.3-ContigExtension/{sample}_mapped.bam",
        CheckCoverM = "Checks/4.3-ContigExtension_{sample}_CoverM.done",
        CheckCE_filtering = "Checks/4.3-ContigExtension_{sample}_filtering.done"
    output:
        CheckCOBRA = "Checks/4.3-ContigExtension_COBRA_{sample}.done",
    threads: 16
    resources:
        mem_mb = 48000,
        cores = 16,
        partition = "high2",
        time = "3-00:00:00"
    params:
        tag = "{sample}",
        output_folder = "4-virus_identification/4.3-ContigExtension/{sample}",
    shell:'''
    micromamba run -n COBRA_env cobra-meta -f {input.Fasta} -q {input.Query} -c {input.CoverageDoc} -m {input.Bam} -o {params.output_folder}_COBRA -a megahit -mink 27 -maxk 127 \
    && \
    touch {output.CheckCOBRA}
    '''

    #### 4.4 - Genomad Rerun ####
rule GeNomadRerun:
    input:
        CheckCobra = "Checks/4.3-ContigExtension_COBRA_{sample}.done",
    output:
        CheckGenomadRerun = "Checks/4.4_GenomadRerun_{sample}.done"
    threads: 16
    resources:
        mem_mb = 64000,
        partition = "high2",
        time = "3-00:00:00"
    params:
        tag = "{sample}_rerun",
        GenomadDB = "/group/jbemersogrp/databases/genomad/genomad_db",
        Contigs = "4-virus_identification/4.3-ContigExtension/{sample}_COBRA/{sample}_renamed_contigs.new.fa"
    shell:'''
        mkdir -p 4-virus_identification/genomad_rerun/{params.tag} && \
        micromamba run -n genomad_env genomad end-to-end --cleanup --composition virome --enable-score-calibration -t {threads} {params.Contigs} \
        4-virus_identification/genomad_rerun/{params.tag} {params.GenomadDB} \
        && \
        touch {output.CheckGenomadRerun}
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
        CheckRenameAssemblyContigs = expand("Checks/3.4-RenameContigs_{sample}.done", sample=rna_viromes),
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