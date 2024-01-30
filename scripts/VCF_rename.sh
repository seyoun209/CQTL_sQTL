#!/bin/bash

ml plink/1.90b3
ml samtools/1.9
ml vcftools/0.1.15

bcftools reheader --sample $2 $1 -o $3
tabix -p vcf $3
bcftools annotate $3 --rename-chrs /work/users/s/e/seyoun/CQTL_sQTL/scripts/chrnm.txt -Oz -o $4
tabix -p vcf $4

