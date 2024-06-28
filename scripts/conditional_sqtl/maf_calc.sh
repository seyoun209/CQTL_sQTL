#!/bin/bash
  
ml plink/1.90b3
ml samtools/1.9
ml vcftools/0.1.15

plink --vcf /work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/06.subset_sigSNps/sigSNPs_all_pbs_vcf.gz \
	--freq \
	--double-id \
	--out /work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/06.subset_sigSNps/sigSNPs_freq

plink --vcf /work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/06.subset_sigSNps/sigSNPs_all_fnf_vcf.gz \
        --freq \
	--double-id \
        --out /work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/06.subset_sigSNps/sigSNPs_freq
