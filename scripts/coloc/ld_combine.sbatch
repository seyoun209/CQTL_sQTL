#!/bin/bash
#SBATCH -J qtltools
#SBATCH -p general
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --mem=30g
#SBATCH -t 3-00:00:00
#SBATCH -e ../../output/logs/combined_ld_qtls_e.%j
#SBATCH -o ../../output/logs/combined_ld_qtls_o.%j

cd /work/users/s/e/seyoun/CQTL_sQTL

ml r/4.3.1



# Run the R script with the appropriate arguments
Rscript /work/users/s/e/seyoun/CQTL_sQTL/scripts/coloc/ld05_combined_sQTL.R "$1" "$2" "$3"
