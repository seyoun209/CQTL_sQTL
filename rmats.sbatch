#!/bin/bash
#SBATCH -J rmats
#SBATCH -p general
#SBATCH -n 12
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -t 3-00:00:00
#SBATCH -e ./rmats_e.%j
#SBATCH -o ./rmats_o.%j

module load rmats-turbo/4.1.1

mkdir -p output/rmats_fnf
mkdir -p output/rmats_oa

run_rmats --b1 output/align/ctl.txt --b2 output/align/fnf.txt \
	--gtf /work/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.gtf \
	-t paired --readLength 100 --nthread 12 \
	--od output/rmats_fnf \
       	--tmp output/rmats_fnf/rmats_ctl_fnf

run_rmats --b1 output/align/ctl.txt --b2 output/align/oa.txt \
	--gtf /work/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.gtf \
	-t paired --readLength 100 --nthread 12 \
	--od output/rmats_oa \
	--tmp output/rmats_oa/rmats_ctl_oa
