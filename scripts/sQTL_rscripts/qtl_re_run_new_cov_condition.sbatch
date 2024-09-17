#!/bin/bash
#SBATCH -J qtl_perm
#SBATCH -p general
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --mem=10g
#SBATCH -t 20:00:00
#SBATCH -e /work/users/s/e/seyoun/CQTL_sQTL/output/logs/qtltools_condition_to_cov_e.%j
#SBATCH -o /work/users/s/e/seyoun/CQTL_sQTL/output/logs/qtltools_condition_to_cov_o.%j


ml qtltools/1.3.1


CHR=$1

# Base directory for output
BASE_DIR="/work/users/s/e/seyoun/CQTL_sQTL/output/02.qtltools_condition_cov"

# Create necessary directories
mkdir -p ${BASE_DIR}/{nominal,perm}

for i in {1..20}
do
    echo "Processing PC${i} for chromosome ${CHR}"

    # Create PC-specific directories
    mkdir -p ${BASE_DIR}/nominal/pc${i}
    mkdir -p ${BASE_DIR}/perm/pc${i}

QTLtools cis --vcf /work/users/s/e/seyoun/CQTL_sQTL/output/geno/merged/merged_pbs_fnf_output.vcf.gz --bed /work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster/qtltools_prep/qtl_condition_to_covariates/ctrvsfnf_qqnorm_chr${CHR}.bed.gz --cov /work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster/qtltools_prep/qtl_condition_to_covariates/covariates_PC${i} --window 100000 --out /work/users/s/e/seyoun/CQTL_sQTL/output/02.qtltools_condition_cov/nominal/pc${i}/chr${CHR}.cis --nominal 1.0 --std-err

QTLtools cis --vcf /work/users/s/e/seyoun/CQTL_sQTL/output/geno/merged/merged_pbs_fnf_output.vcf.gz --bed /work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster/qtltools_prep/qtl_condition_to_covariates/ctrvsfnf_qqnorm_chr${CHR}.bed.gz --cov /work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster/qtltools_prep/qtl_condition_to_covariates/covariates_PC${i}  --window 100000 --out /work/users/s/e/seyoun/CQTL_sQTL/output/02.qtltools_condition_cov/perm/pc${i}/chr${CHR}.perm --permute 1000 --std-err

done
