#!/bin/bash

# define the directory that contains the files
directory=$1

# define the output file
output_file=$2

# remove the output file if it already exists
if [ -f "$output_file" ]; then
    rm "$output_file"
fi

# iterate over the files in the directory
for file in "$directory"/*_R1_001.fastq.gz
do
    # extract the filename without the extension
    filename=$(basename -- "$file")
    filename="${filename%_R1_001.fastq.gz}"

    # append the filename to the output file
    echo "$filename" >> "$output_file"
done
