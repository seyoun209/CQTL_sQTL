#!/bin/bash
#SBATCH -J samtools
#SBATCH -p general
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -t 10:00:00
#SBATCH -e /work/users/s/e/seyoun/CQTL_sQTL/output/logs/samtools_snrnp70_subset_e.%j
#SBATCH -o /work/users/s/e/seyoun/CQTL_sQTL/output/logs/samtools_snrnp70_subset_o.%j


module load samtools/1.9

#cd /work/users/s/e/seyoun/crispr/02.test_seq

samtools view -bh $1 chr19:49085419-49108605 > $2
samtools view $2 > $3
samtools index $2 $4

