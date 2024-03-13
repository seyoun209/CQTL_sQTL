#!/bin/bash
#SBATCH -J mergesignal
#SBATCH -p general
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=10g
#SBATCH -t 11-00:00:00
#SBATCH -e ./ctl_merge_e.%j
#SBATCH -o ./ctl_merge_o.%j


module load samtools/1.9
module load deeptools/3.5.4

mkdir -p output/signals/05.mergesignal

bamCoverage -b output/signals/03.mergeAlign/CTL_sorted.bam -o output/signals/05.mergesignal/CTL.bw > output/logs/mergeSignal_CTL.err 2>&1

#samtools index output/signals/03.mergeAlign/CTL_sorted.bam >> output/logs/mergeAlign_CTL.err 2>&1

#samtools flagstat output/signals/04.mergeAlign_wasp/CTL_sorted.bam > output/signals/04.mergeAlign_wasp/CTL_stats.txt >> output/logs/mergeAlign_wasp_CTL.err 2>&1
#samtools index output/signals/04.mergeAlign_wasp/CTL_sorted.bam >> output/logs/mergeAlign_wasp_CTL.err 2>&1

