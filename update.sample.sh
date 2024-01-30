#!/bin/bash

#trim_directory='/work/users/s/e/seyoun/CQTL_sQTL/output/trim'
fastq_directory="/work/users/s/e/seyoun/CQTL_sQTL/output/fastq"
awk -F'\t' -v fastq_directory="$fastq_directory" 'BEGIN {OFS=FS} {
    trim1_file = fastq_directory "/" $1 "_" $4 "_" $6 "_" $11 "_R1.fastq.gz"
    trim2_file = fastq_directory "/" $1 "_" $4 "_" $6 "_" $11 "_R2.fastq.gz"

    # Update the line with trim1 and trim2 columns
    $12 = trim1_file
    $13 = trim2_file

    # Print the updated line
    print
}' samplesheet.txt | awk -F'\t' '!seen[$1,$4,$6,$12,$13]++' > temp.txt


awk -F'\t' 'BEGIN {OFS=FS} NR==1 {$12="cat1"; $13="cat2"} 1' temp.txt > updated_samplesheet.txt

rm temp.txt


