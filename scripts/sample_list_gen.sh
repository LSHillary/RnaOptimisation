#!/bin/bash

# define the directory that contains the files
directory=$1

# define the output file
output_file="samples.txt"

# remove the output file if it already exists
if [ -f "$output_file" ]; then
    rm "$output_file"
fi

# iterate over the files in the directory
for file in "$directory"/*_raw_R1.fq.gz
do
    # extract the filename without the extension
    filename=$(basename -- "$file")
    filename="${filename%_raw_R1.fq.gz}"

    # append the filename to the output file
    echo "$filename" >> "$output_file"
done
