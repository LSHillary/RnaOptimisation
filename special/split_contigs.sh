#!/bin/bash

# Create a new directory to store all the contig FASTA files
mkdir -p contigs_combined

# Initialize a counter for contig renaming
counter=1

# Loop through each .fna file in the current directory
for fna_file in *.fna; do
    # Generate the new filename with the counter
    new_filename="contigs_combined/contig_${counter}.fna"

    # Rename the .fna file to the new filename
    mv "$fna_file" "$new_filename"

    # Increment the counter for the next contig
    ((counter++))

    # Generate a file-specific prefix for contig headers
    prefix="${new_filename%.fna}"

    # Split the .fna file into individual contigs using seqtk
    seqtk seq -A "$new_filename" |
    awk -v filename="$prefix" '/^>/{out=filename"_"substr($0,2)".fasta"; print > out; next}{print >> out}'
done
