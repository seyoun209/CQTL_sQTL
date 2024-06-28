#!/bin/bash
#BATCH -J plink
#SBATCH -p general
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=20g
#SBATCH -t 2-00:00:00
#SBATCH -e ../../output/logs/ld_calc_1000g_e.%j
#SBATCH -o ../../output/logs/ld_calc_1000g_o.%j

ml plink/1.90b3
ml qtltools/1.3.1
ml samtools/1.20


r2=0
window_kb=1000
ld_window=200000


ref_our=/work/users/s/e/seyoun/CQTL_sQTL/output/geno/vcf_bed/$1


mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld/$1
mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld05/$1
mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld_1000g_EUR/$1
mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld05_1000g_EUR/$1
mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld_1000g_ALL/$1
mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld05_1000g_ALL/$1


#ref_fnf="/work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz"

plink --bfile ${ref_our} \
      --ld-snp $2 \
      --ld-window ${ld_window} \
      --ld-window-kb ${window_kb} \
      --ld-window-r2 0 \
      --r2 \
      --threads 1 \
      --memory 4098 \
      --out /work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld/$1/$2

plink --bfile ${ref_our} \
      --ld-snp $2 \
      --ld-window ${ld_window} \
      --ld-window-kb ${window_kb} \
      --ld-window-r2 0.5 \
      --r2 \
      --threads 1 \
      --memory 4098 \
      --out /work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld05/$1/$2

plink --bfile $3 \
      --ld-snp $2 \
      --ld-window ${ld_window} \
      --ld-window-kb ${window_kb} \
      --ld-window-r2 0 \
      --r2 \
      --threads 1 \
      --memory 4098 \
      --out /work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld_1000g_EUR/$1/$2

plink --bfile $3 \
      --ld-snp $2 \
      --ld-window ${ld_window} \
      --ld-window-kb ${window_kb} \
      --ld-window-r2 0.5 \
      --r2 \
      --threads 1 \
      --memory 4098 \
      --out /work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld05_1000g_EUR/$1/$2


plink --bfile $4 \
      --ld-snp $2 \
      --ld-window ${ld_window} \
      --ld-window-kb ${window_kb} \
      --ld-window-r2 0 \
      --r2 \
      --threads 1 \
      --memory 4098 \
      --out /work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld_1000g_ALL/$1/$2

plink --bfile $4 \
      --ld-snp $2 \
      --ld-window ${ld_window} \
      --ld-window-kb ${window_kb} \
      --ld-window-r2 0.5 \
      --r2 \
      --threads 1 \
      --memory 4098 \
      --out /work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld05_1000g_ALL/$1/$2


