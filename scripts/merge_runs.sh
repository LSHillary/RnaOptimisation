#!/bin/bash

# Check if the right number of arguments are provided
if [ $# -ne 3 ]; then
  echo "Usage: $0 <dir1> <dir2> <merged_dir>"
  exit 1
fi

# Directories
dir1=$1
dir2=$2
merged_dir=$3

# Create the merged directory if it doesn't exist
mkdir -p $merged_dir

# Get all file names in the first directory
for file in "$dir1"/*; do
    # Extract the base name of the file
    base_name=$(basename "$file")

    # If the file exists in both directories
    if [ -f "$dir2/$base_name" ]; then
        # Merge the files
        cat "$dir1/$base_name" "$dir2/$base_name" > "$merged_dir/$base_name"
        echo "Merged $base_name"
    fi
done
