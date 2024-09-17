#/bin/bash
#SBATCH -J qtl_perm
#SBATCH -p general
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --mem=10g
#SBATCH -t 10:00:00
#SBATCH -e /work/users/s/e/seyoun/CQTL_sQTL/output/logs/qtltools_nominal_1mb_e.%j
#SBATCH -o /work/users/s/e/seyoun/CQTL_sQTL/output/logs/qtltools_nominal_1mb_o.%j


ml qtltools/1.3.1

#for i in {1..20}
#do
#echo ${i}

mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/nominal_1mb/nominal_pbs/pc5
mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/nominal_1mb/nominal_fnf/pc4


