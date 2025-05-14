#!/bin/bash
#SBATCH -J runFDR_1
#SBATCH -p general
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=10g
#SBATCH -t 2-00:00:00
#SBATCH -e ../../output/logs/qtl_sig_e.%j
#SBATCH -o ../../output/logs/qtl_sig_o.%j

ml r/4.4.0

Rscript runFDR_cis.R /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/perm_pbs/pc5_allchr.pbs.perm 1 /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/01.significant/pbs_1_pc5

Rscript runFDR_cis.R /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/perm_fnf/pc4_allchr.fnf.perm 1 /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/01.significant/fnf_1_pc4


