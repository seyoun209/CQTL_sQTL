#!/bin/bash
#SBATCH -J rsid
#SBATCH -p general
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -t 10:00:00
#SBATCH -e ../../output/logs/rsid_chr$1.%j.err
#SBATCH -o ../../output/logs/rsid_chr$1.%j.out

ml samtools/1.20

dir="/work/users/s/e/seyoun/CQTL_sQTL/output/geno/"
chr=$1

echo "Processing chromosome $1"
output_file="${dir}rsid/chr$1_snps_info_with_rsid.bed"
temp_dbsnp="${dir}chr$1_temp_dbsnp.txt"
temp_input="${dir}chr$1_temp_input.txt"
dbsnp_vcf="/users/s/e/seyoun/09.cqtl_hic/oa_tissue/RNAPIPE-Splice/output/Ref/dbSNP/dbSNP155.GRCh38.p13_chr$1.vcf.gz"

mkdir -p "${dir}rsid"

awk -v chr="chr$1" '$1 == chr' "${dir}cqtl_snps_info.txt" > "$temp_input"
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "$dbsnp_vcf" > "$temp_dbsnp"

# Process input file and add dbSNP IDs
awk 'BEGIN {OFS="\t"}
     NR==FNR {
       split($5, alts, ",");
       for (i in alts) {
         key = $1 ":" $2 ":" $4 ":" alts[i];
         a[key] = $3;
         key = $1 ":" $2 ":" alts[i] ":" $4;
         a[key] = $3;
       }
       next
     }
     {
       key = $1 ":" $2 ":" $4 ":" $5;
       if (key in a) {
         print $0, a[key];
       } else {
         print $0, ".";
       }
     }' "$temp_dbsnp" "$temp_input" > "$output_file"

# Clean up temporary files
#rm "$temp_dbsnp" "$temp_input"

echo "Finished processing chromosome $1"
