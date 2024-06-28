#!/bin/bash
#SBATCH -J mergesignal
#SBATCH -p general
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=60g
#SBATCH -t 2-00:00:00
#SBATCH -e ../../output/logs/response_qtl_e.%j
#SBATCH -o ../../output/logs/response_qtl_o.%j


ml r/4.3.1

Rscript /work/users/s/e/seyoun/CQTL_sQTL/scripts/sQTL_rscripts/response_qtl_v2.R
