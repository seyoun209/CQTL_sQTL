#!/bin/bash
#SBATCH -J qtltools
#SBATCH -p general
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=20g
#SBATCH -t 3-00:00:00
#SBATCH -e /work/users/s/e/seyoun/CQTL_sQTL/output/logs/editing_ld0_sqtl_e.%j
#SBATCH -o /work/users/s/e/seyoun/CQTL_sQTL/output/logs/editing_ld0_sqtl_o.%j

ml r/4.3.1

Rscript /work/users/s/e/seyoun/CQTL_sQTL/scripts/coloc/ld05_combined_sQTL.R 

