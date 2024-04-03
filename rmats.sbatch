#!/bin/bash
#SBATCH -J rmats
#SBATCH -p general
#SBATCH -n 12
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -t 3-00:00:00
#SBATCH -e ./output/logs/rmats_e.%j
#SBATCH -o ./output/logs/rmats_o.%j

module load rmats-turbo/4.1.1

ml r/4.3.1

#Rscript /work/users/s/e/seyoun/CQTL_sQTL/scripts/Differential_splicing/snrnp70_het_kd_rmats_sigAS.R 

mkdir -p output/rmats_fnf_101
mkdir -p output/rmats_oa_101

run_rmats --b1 output/align/ctl.txt --b2 output/align/fnf.txt \
	--gtf /work/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.gtf \
	-t paired --readLength 100 --nthread 12 \
	--od output/rmats_fnf_101 \
       	--tmp output/rmats_fnf_101/rmats_ctl_fnf

run_rmats --b1 output/align/ctl.txt --b2 output/align/oa.txt \
	--gtf /work/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.gtf \
	-t paired --readLength 100 --nthread 12 
	--od output/rmats_oa_101 \
	--tmp output/rmats_oa_101/rmats_ctl_oa
Rscript /work/users/s/e/seyoun/CQTL_sQTL/scripts/Differential_splicing/snrnp70_het_kd_rmats_sigAS.R


