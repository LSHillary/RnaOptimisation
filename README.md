# ViromeDataProcessing
Central codebase for virome data processing pipeline

## Step 0 - Set up and fetch raw data

### Installing Conda

Personally, I'm going to try micromamba:

https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html

`curl micro.mamba.pm/install.sh | bash`

### Setting up the ViromeProcessing conda environment

`micromamba create -n ViromeDataProcessing -f ViromeDataProcessing.yml`

### Non-conda resources

Follow these instructions to install bbtools
`https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/`

## Step 0 - Get data and verify

Run the following command:

`conda activate ViromeDataProcessing`

### Download raw data
For people wanting to reproduce this analysis post publication, download raw reads from the SRA

For me, I copied the data from the UC Davis Genome Center server into

### Run 0-preprocessing.smk
This has the option of merging multiple runs and runs FastQC and MultiQC 
