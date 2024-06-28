#!/bin/bash
#SBATCH -J qtltools
#SBATCH -p general
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --mem=10g
#SBATCH -t 10:00:00
#SBATCH -e output/logs/cond_re_run_e.%j
#SBATCH -o output/logs/cond_re_run_o.%j

cd /work/users/s/e/seyoun/CQTL_sQTL

export PATH=/nas/longleaf/rhel8/apps/qtltools/1.3.1/bin:/nas/longleaf/rhel8/apps/samtools/1.20/bin:/nas/longleaf/rhel8/apps/samtools/1.20/htslib/bin:/nas/longleaf/home/seyoun/miniconda3/bin:/nas/longleaf/apps/tabix/0.2.6/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/puppetlabs/bin:/nas/longleaf/home/seyoun/.local/bin:/nas/longleaf/home/seyoun/bin:$PATH


echo $PATH
which QTLtools

vcf="$1"
bed="$2"
norm_out="$3"
perm_out="$4"
pheno="$5"
cov="$6"
pheid="$7"

echo "VCF file: $vcf"
echo "Norm Dir: $norm_out"
echo "Perm Dir: $perm_out"
echo "Beddir: $bed"
echo "Cov Dir: $cov"
echo "Pheno_txt: $pheno"
echo "Phenotype ID: $pheid"

QTLtools cis --vcf $vcf --bed $bed --cov $cov --window 100000 --include-phenotypes $pheno --out $norm_out --nominal 1.0 --std-err

QTLtools cis --vcf $vcf --bed $bed --cov $cov --window 100000 --include-phenotypes $pheno --out $perm_out --permute 1000 

