#!/bin/bash

ml samtools/1.9

bcftools_output=$(bcftools query -l $1)

# Loop through each line in the bcftools output
while read -r line; do
    # Extract the string after "AM" in the current line
    sample_name=$(echo "$line" | grep -oP 'AM[0-9]{4}')

    # Find the corresponding line in tissue_meta.txt
    matched_sample=$(grep -F "$sample_name" $2 |awk '{print $1}')
    # Print the sample name and the matched sample name 

    printf "%-25s %s\n" "$line" "$matched_sample"
done <<< "$bcftools_output" 


