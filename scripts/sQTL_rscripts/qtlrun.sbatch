#!/bin/bash
#SBATCH -J qtl_perm
#SBATCH -p general
#SBATCH -n 12
#SBATCH -N 1
#SBATCH --mem=20g
#SBATCH -t 20:00:00
#SBATCH -e ./tmp/qtltools_e.%j
#SBATCH -o ./tmp/qtltools_o.%j


ml qtltools/1.3.1

for i in {1..20}
do
echo ${i}

QTLtools cis --vcf /work/users/s/e/seyoun/07.Doug_thesis/09.cqtl_hic/oa_tissue/geno/pbs_geno/04.final/pbs_snpfiltered.fi_withchr_noindel.vcf.gz --bed /work/users/s/e/seyoun/07.Doug_thesis/09.cqtl_hic/oa_tissue/RNAPIPE-Splice/output/clu_fnf/ctrvsfnf_qqnorm_chr$1.bed.gz --cov ./clu_fnf/covariates_final/nobatch_cov/nobatch_covariates_PC${i} --window 100000 --out 04.qtl_noindelbatch/nominal_pbs/pc${i}/chr$1.pbs.cis --nominal 1.0

QTLtools cis --vcf /work/users/s/e/seyoun/07.Doug_thesis/09.cqtl_hic/oa_tissue/geno/fnf_geno/04.final/fnf_snpfiltered.fi_withchr_noindel.vcf.gz --bed /work/users/s/e/seyoun/07.Doug_thesis/09.cqtl_hic/oa_tissue/RNAPIPE-Splice/output/clu_fnf/ctrvsfnf_qqnorm_chr$1.bed.gz --cov ./clu_fnf/covariates_final/nobatch_cov/nobatch_covariates_PC${i} --window 100000  --out 04.qtl_noindelbatch/nominal_fnf/pc${i}/chr$1.fnf.cis --nominal 1.0


QTLtools cis --vcf /work/users/s/e/seyoun/07.Doug_thesis/09.cqtl_hic/oa_tissue/geno/pbs_geno/04.final/pbs_snpfiltered.fi_withchr_noindel.vcf.gz --bed /work/users/s/e/seyoun/07.Doug_thesis/09.cqtl_hic/oa_tissue/RNAPIPE-Splice/output/clu_fnf/ctrvsfnf_qqnorm_chr$1.bed.gz --cov ./clu_fnf/covariates_final/nobatch_cov/nobatch_covariates_PC${i} --window 100000 --out 04.qtl_noindelbatch/perm_pbs/pc${i}/chr$1.pbs.perm --permute 1000

QTLtools cis --vcf /work/users/s/e/seyoun/07.Doug_thesis/09.cqtl_hic/oa_tissue/geno/fnf_geno/04.final/fnf_snpfiltered.fi_withchr_noindel.vcf.gz --bed /work/users/s/e/seyoun/07.Doug_thesis/09.cqtl_hic/oa_tissue/RNAPIPE-Splice/output/clu_fnf/ctrvsfnf_qqnorm_chr$1.bed.gz --cov ./clu_fnf/covariates_final/nobatch_cov/nobatch_covariates_PC${i} --window 100000  --out 04.qtl_noindelbatch/perm_fnf/pc${i}/chr$1.fnf.perm --permute 1000


done
