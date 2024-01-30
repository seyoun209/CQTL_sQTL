#!/bin/bash

ml samtools/1.9

bcftools_output=$(bcftools query -l $1)

# Loop through each line in the bcftools output
#while read -r line; do
#    # Extract the string after "AM" in the current line
#    sample_name=$(echo "$line" | grep -oP 'AM[0-9]{4}')
#
#    # Find the corresponding line in tissue_meta.txt
#    matched_sample=$(grep -F "$sample_name" $2 |awk '{print $1}')
#    # Print the sample name and the matched sample name 
#
#    printf "%-25s %s\n" "$line" "$matched_sample"
#done <<< "$bcftools_output" 


while read -r line; do
    # Check if the line contains "AM" (PBS/FNF) or "OA" (OA)
    if [[ "$line" =~ AM[0-9]{4} ]]; then
        # Extract the string after "AM" in the current line
        sample_name=$(echo "$line" | grep -oP 'AM[0-9]{4}')
    elif [[ "$line" =~ OA[0-9]{3} ]]; then
        # Extract the string after "OA" in the current line
        sample_name=$(echo "$line" | grep -oP 'OA[0-9]{3}')
    else
        # Handle other cases or skip the line if needed
        continue
    fi

    # Find the corresponding line in tissue_meta.txt
    matched_sample=$(awk -v sample="$sample_name" '$1 ~ sample {print $1}' "$2")

    # Print the sample name and the matched sample name
    printf "%-25s %s\n" "$line" "$matched_sample"

done <<< "$bcftools_output" > $3


