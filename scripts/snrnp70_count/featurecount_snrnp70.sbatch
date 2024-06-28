#!/bin/bash
#SBATCH -J featurecount
#SBATCH -p general
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem=40g
#SBATCH -t 1-00:00:00
#SBATCH -e ../../output/logs/featurecount_snrnp70_e.%j
#SBATCH -o ../../outout/logs/featurecount_snrnp70_o.%j


ml subread/2.0.6

dir_bam="/work/users/s/e/seyoun/CQTL_sQTL/output/snrnp70_subset_bam"
featureCounts -a /work/users/s/e/seyoun/Ref/genome/snrnp70_flatten_subset.gtf \
        -o /work/users/s/e/seyoun/CQTL_sQTL/output/featurecounts_snrnp70/snrnp70counts.txt \
        -R BAM \
        -T 8 \
        -g exon_id \
        -p \
        -O \
        -F SAF \
        -t exon \
        ${dir_bam}/*_sorted.bam_snrnp70.bam > ../../output/logs/featurecount.err 2>&1


