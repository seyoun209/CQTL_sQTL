#!/bin/bash

# Define the directory containing the sample files
sample_directory="/work/users/s/e/seyoun/CQTL_sQTL/output/fastq"

# Get the sample names from the files in the directory
sample_names=("$sample_directory"/*_R1.fastq.gz)
sample_names=${sample_names%%_R1.fastq.gz}

# Define mapping for categories
declare -A category_mapping=(
    ["CTL"]="CTL"
    ["FNF"]="FNF"
    ["OA"]="OA"
)

# Initialize arrays for each category
declare -A category_samples

# Categorize samples
for sample_path in "${sample_names[@]}"; do
    # Extract the sample name from the file path
    sample_name=$(basename "$sample_path")
    sample_name=${sample_name%%_R1.fastq.gz}

    for category in "${!category_mapping[@]}"; do
        pattern="${category_mapping[$category]}"
        if [[ $sample_name == *$pattern* ]]; then
            category_samples["$category"]+="$sample_name $category"$'\n'
        fi
    done
done

# Save categorized samples into separate files
for category in "${!category_samples[@]}"; do
    output_file="/work/users/s/e/seyoun/CQTL_sQTL/output/${category}_samples.txt"
    printf "%s\n" "${category_samples[$category]}" > "$output_file"
done

