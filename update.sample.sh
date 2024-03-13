#!/bin/bash

#trim_directory='/work/users/s/e/seyoun/CQTL_sQTL/output/trim'
fastq_directory="/work/users/s/e/seyoun/CQTL_sQTL/output/fastq"
align_dir='/work/users/s/e/seyoun/CQTL_sQTL/output/align'

awk -F'\t' -v dir_test="$align_dir" 'BEGIN {OFS=FS} {
    file1 = dir_test "/" $1 "_" $4 "_" $6 "_" $11 "_sorted.bam"
    #file2 = dir_test "/" $1 "_" $4 "_" $6 "_" $11 "_R2.fastq.gz"

    # Update the line with trim1 and trim2 columns
    $12 = file1
    #$13 = file2

    # Print the updated line
    print
}' samplesheet.txt | awk -F'\t' '!seen[$1,$4,$6,$12]++' > temp.txt

awk -F'\t' 'BEGIN {OFS=FS} NR==1 {$12="alignbam"} 1' temp.txt > temp_samplesheet.txt
#awk -F'\t' 'BEGIN {OFS=FS} NR==1 {$12="cat1"; $13="cat2"} 1' temp.txt > updated_samplesheet.txt

rm temp.txt


align_wasp_dir='/work/users/s/e/seyoun/CQTL_sQTL/output/align_wasp'

awk -F'\t' -v dir_test="$align_wasp_dir" 'BEGIN {OFS=FS} {
    file1 = dir_test "/" $1 "_" $4 "_" $6 "_" $11 ".Aligned.sortedByCoord.WASP.bam"
    

    # Update the line with trim1 and trim2 columns
    $13 = file1
   

    # Print the updated line
    print
}' temp_samplesheet.txt | awk -F'\t' '!seen[$1,$4,$6,$12,$13]++' > temp.txt

awk -F'\t' 'BEGIN {OFS=FS} NR==1 {$13="alignbam_wasp"} 1' temp.txt > aligned_samplesheet.txt
#awk -F'\t' 'BEGIN {OFS=FS} NR==1 {$12="cat1"; $13="cat2"} 1' temp.txt > updated_samplesheet.txt

rm temp.txt


