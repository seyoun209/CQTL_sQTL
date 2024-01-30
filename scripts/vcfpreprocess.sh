#!/bin/bash

ml plink/1.90b3
ml samtools/1.9
ml vcftools/0.1.15

bcftools reheader --sample $2 $1 -o $3 
tabix -p vcf $3
bcftools annotate $3 --rename-chrs /work/users/s/e/seyoun/CQTL_sQTL/scripts/chrnm.txt -Oz -o $4
tabix -p vcf $4

#make bed bim fam files
mkdir -p $5/01.bed
plink --vcf $4 \
        --make-bed \
        --double-id \
        --out $5/01.bed/cqtl

# calculate allele frequency and filter minor allel 10 and minHets: 5
mkdir -p 0$5/2.freq_files

plink --bfile $5/01.bed/cqtl\
        --freqx \
        --out $5/02.freq_files/freq

awk 'NR==1 || ($6 >= 5  && ($5 >= 10 || ($5 + $7) >= 10))' $5/02.freq_files/freq.frqx > $5/02.freq_files/snp_list.txt;
sed -i '1d' $5/02.freq_files/snp_list.txt;


# filter out snp list
mkdir -p $5/03.recoded
plink --bfile $5/01.bed/cqtl \
        --extract $5/02.freq_files/snp_list.txt \
        --make-bed \
        --recode 12 \
        --output-missing-genotype 0 \
        --transpose \
        --out $5/03.recoded/recoded

# change back to vcf file
mkdir -p $5/04.final
plink --bfile $5/03.recoded/recoded_fnf \
        --recode vcf \
        --out $5/04.final/finalFiltered

bgzip $5/04.final/finalFiltered.vcf

bcftools query -l  $5/04.final/finalFiltered_fnf.vcf.gz | awk -F_ '{print $0"\t"$1"_"$2"_"$3}' > $5/04.final/sample_nm.txt
bcftools reheader --sample $5/04.final/sample_nm.txt $5/04.final/finalFiltered_fnf.vcf.gz -o $5/04.final/snpfiltered_renamed.vcf.gz

rm -rf $5/04.final/finalFiltered_fnf.vcf.gz

tabix -p vcf $5/04.final/fnf_snpfiltered.fi.vcf.gz

vcftools --vcf $5/04.final/snpfiltered_renamed.vcf.gz --remove-indels --recode --recode-INFO-all --out $5/04.final/snpfiltered_wCHR_no_indels.final.vcf

bgzip $5/04.final/snpfiltered_wCHR_no_indels.final.vcf 
tabix -p vcf $5/04.final/snpfiltered_wCHR_no_indels.final.vcf.gz
bcftools stats $5/04.final/snpfiltered_wCHR_no_indels.final.vcf.gz > $5/04.final/snpfiltered_wCHR_no_indels.final_stats.txt


#creating PCA

plink --vcf $5/04.final/snpfiltered_wCHR_no_indels.final.vcf --pca  --out $5/05.pca/CQTL --double-id











