#!/bin/bash
#SBATCH -J deseq
#SBATCH -p general
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=30g
#SBATCH -t 3-00:00:00
#SBATCH -e ../output/logs/deseq_e.%j
#SBATCH -o ../output/logs/deseq_o.%j


ml r/4.3.1


Rscript diff_transcript_dds.R
