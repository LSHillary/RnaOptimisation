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

DnaRuns = extract_items(data, ['DnaRuns'])

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
        #expand("Checks/4.0-GenomadPrefilter_{sample}.done", sample = dna_viromes),
        CheckGenomadRerun = expand("Checks/4.4_GenomadRerun_{sample}.done", sample = dna_viromes),

# Run GeNomad on coassembled contigs
rule GeNomadIndividual:
    input:
        ContigsIn = "3-Assembly/Contigs/{sample}_renamed_contigs.fna",
    output:
        CheckGenomad = "Checks/4.1-Genomad_{sample}.done"
    threads: 12
    params:
        tag = "{sample}",
        GenomadDB = "/group/jbemersogrp/databases/genomad/genomad_db",
        GenomadFolder = "4-virus_identification/genomad_all/{sample}"
    resources:
        mem_mb = 65536,
        partition = "high2",
        time = "3-00:00:00"
    shell:'''
        mkdir -p {params.GenomadFolder} && \
        micromamba run -n genomad_env genomad end-to-end --cleanup --composition virome --enable-score-calibration -t {threads} {input.ContigsIn} \
        {params.GenomadFolder} {params.GenomadDB} && \
        touch {output.CheckGenomad}
        '''

#### 4.2 - Contig Extension ####
# 4.2.1 - Mapping
rule CE_Mapping:
    input:
        FwdReadsRun1="2-QC/2.4-ErrorCorrection/{sample}_run1_EC_R1.fq.gz",
        RevReadsRun1="2-QC/2.4-ErrorCorrection/{sample}_run1_EC_R2.fq.gz",
        FwdReadsRun2="2-QC/2.4-ErrorCorrection/{sample}_run2_EC_R1.fq.gz",
        RevReadsRun2="2-QC/2.4-ErrorCorrection/{sample}_run2_EC_R2.fq.gz",
        FwdReadsRun3="2-QC/2.4-ErrorCorrection/{sample}_run3_EC_R1.fq.gz",
        RevReadsRun3="2-QC/2.4-ErrorCorrection/{sample}_run3_EC_R2.fq.gz",
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
        cat {input.FwdReadsRun1} {input.FwdReadsRun2} {input.FwdReadsRun3} > {output.FwdReadsAll} && \
        cat {input.RevReadsRun1} {input.RevReadsRun2} {input.RevReadsRun3} > {output.RevReadsAll} && \
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
        GenomadFasta = "4-virus_identification/genomad_all/{sample}/{sample}_renamed_contigs_summary/{sample}_renamed_contigs_virus.fna"
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

rule FilterGenomadRerun:
    input:
        CheckGenomadRerun = "Checks/4.3-ContigExtension_COBRA_{sample}.done"
    output:
        FilteredGenomadResults = "4-virus_identification/genomad_rerun/{sample}/{sample}_FilteredGenomadResults.tsv",
        CheckFilteredGenomadResults = "Checks/4.3-FilterGenomadResults_{sample}.done",
        FilteredContigs = "4-virus_identification/genomad_rerun/{sample}_FilteredContigs.fna"
    params:
        tag = "{sample}",
        Contigs = "4-virus_identification/genomad_rerun/{sample}/{sample}_summary/{sample}_virus.fna",
        GenomadResults = "4-virus_identification/genomad_rerun/{sample}/{sample}_summary/{sample}_virus_summary.tsv"
    shell:'''
        awk -F'\t' 'NR==1 || $2 > 10000 || $11 ~ /Monodnaviria/' {params.GenomadResults} | \
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