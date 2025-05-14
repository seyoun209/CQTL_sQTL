#!/bin/bash
  
ml plink/2.00a-20220129
ml samtools/1.21
ml bcftools/1.12  # if needed

mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/07.subset_sigSNps_eqtl
mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/07.subset_sigSNps_eqtl

# Define input VCF file paths for each condition
pbs_vcf="/work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz"
fnf_vcf="/work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz"
extract_file="/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/e_sqtl_variant_ids.txt"

# Use an associative array for convenience
declare -A vcf_paths
vcf_paths["pbs"]=$pbs_vcf
vcf_paths["fnf"]=$fnf_vcf


for CONDITION in pbs fnf; do
  echo "Processing condition: ${CONDITION}"
  
  input_vcf=${vcf_paths[${CONDITION}]}
  if [ "${CONDITION}" == "pbs" ]; then
    out_dir="/work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/07.subset_sigSNps_eqtl"
  else
    out_dir="/work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/07.subset_sigSNps_eqtl"
  fi

  # Run plink2 using --extract to subset variants using the provided variant list,
  # then recode to A-transpose format with --const-fid 0 to force FID=0.
  plink2 --vcf "${input_vcf}" \
         --extract "${extract_file}" \
         --const-fid 0 \
         --recode A-transpose \
         --keep-allele-order \
         --out "${out_dir}/recodeA_${CONDITION}"
  
  # Remove the unwanted "0_" prefix from sample names in the header
  sed -i '1s/\(^\|\t\)0_/\1/g' "${out_dir}/recodeA_${CONDITION}.traw"
  
  echo "Finished processing ${CONDITION}; output written to ${out_dir}"
done