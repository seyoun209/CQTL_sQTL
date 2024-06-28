#!/bin/bash

ml gatk/4.5.0.0
ml plink/1.90b3
ml samtools/1.9
ml vcftools/0.1.15

mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/06.subset_sigSNps
mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/06.subset_sigSNps

pbs_dir="/work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/06.subset_sigSNps"
fnf_dir="/work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/06.subset_sigSNps"

gatk SelectVariants -V /work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz \
	--keep-ids /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/sigSNPs_all.list \
	-O ${pbs_dir}/sigSNPs_all_pbs_vcf.gz > /work/users/s/e/seyoun/CQTL_sQTL/output/logs/gatk_pbs.err 2>&1

gatk SelectVariants -V /work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz \
        --keep-ids /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/sigSNPs_all.list \
        -O ${fnf_dir}/sigSNPs_all_fnf_vcf.gz > /work/users/s/e/seyoun/CQTL_sQTL/output/logs/gatk_fnf.err 2>&1


tabix -p vcf ${pbs_dir}/sigSNPs_all_pbs_vcf.gz
tabix -p vcf ${fnf_dir}/sigSNPs_all_fnf_vcf.gz



plink --vcf ${pbs_dir}/sigSNPs_all_pbs_vcf.gz \
	--double-id \
	--recode A \
	--out ${pbs_dir}/recodeA_pbs

plink --vcf ${fnf_dir}/sigSNPs_all_fnf_vcf.gz \
        --recode A \
	--double-id \
	--out ${fnf_dir}/recodeA_fnf






