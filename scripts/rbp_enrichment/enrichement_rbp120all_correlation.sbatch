#!/bin/bash
#SBATCH -J enrichement
#SBATCH -p general
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=60g
#SBATCH -t 2-00:00:00
#SBATCH -e /work/users/s/e/seyoun/CQTL_sQTL/output/logs/rbp_correlation_e.%j
#SBATCH -o /work/users/s/e/seyoun/CQTL_sQTL/output/logs/rbp_correlation_o.%j


ml r/4.3.1

Rscript /work/users/s/e/seyoun/CQTL_sQTL/scripts/enrichment/run_correlation_all120rbps.R


