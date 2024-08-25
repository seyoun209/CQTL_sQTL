#!/bin/bash
#SBATCH -J plotting
#SBATCH -p general
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=60g
#SBATCH -t 2-00:00:00
#SBATCH -e ../../output/logs/locuas_box_qtl_e.%j
#SBATCH -o ../../output/logs/locus_box_qtl_o.%j


ml r/4.3.1

Rscript shared_locus_boxPlot.R
