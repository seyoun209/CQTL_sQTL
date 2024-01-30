#!/usr/bin/env python3

import pandas as pd
import glob
import shutil
import os
from utils.namer import namer


## Read in samplesheet
samples = pd.read_csv(config["updatedsamplesheet"], sep='\t')

## Convert all columns to strings
samples = samples.astype(str)

samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)
read1_trim = samples.groupby('mn')['cat1'].apply(list).to_dict()
read2_trim = samples.groupby('mn')['cat2'].apply(list).to_dict()

rule all:
    input:
        expand("output/align_wasp/{sample}_Aligned.sortedByCoord.out.bam", sample=read1_trim.keys()),
        expand("output/align_wasp/{sample}_Log.out", sample=read1_trim.keys()),
        expand("output/align_wasp/{sample}_Aligned.toTranscriptome.out.bam", sample=read1_trim.keys()),
        expand("output/align_wasp/{sample}_Log.progress.out", sample=read1_trim.keys()),
        expand("output/align_wasp/{sample}_STARtmp", sample=read1_trim.keys()),
        expand("output/align_wasp/{sample}_Log.final.out", sample=read1_trim.keys()),
        expand("output/align_wasp/{sample}_SJ.out.tab", sample=read1_trim.keys()),
        expand('output/align_wasp/{sample}_Aligned.sortedByCoord.out.bam.bai',sample=read1_trim.keys())

def generate_vcf_path(sample):
    print(sample)
    if "CTL" in sample:
        print("Matching CTL")
        return "output/geno/pbs_geno/pbs.rename_wCHR.vcf"
    elif "FNF" in sample:
        print("Matching FNF")
        return "output/geno/fnf_geno/fnf.rename_wCHR.vcf"
    elif "OA" in sample:
        print("Matching OA")
        return "output/geno/oa_geno/oa.rename_wCHR.vcf"

rule align:
    input:
        R1 = lambda wildcards: read1_trim.get(wildcards.sample),
        R2 = lambda wildcards: read1_trim.get(wildcards.sample),
        
    output:
        bam = "output/align_wasp/{sample}_Aligned.sortedByCoord.out.bam",
        log_out = "output/align_wasp/{sample}_Log.out",
        transcriptome_bam = "output/align_wasp/{sample}_Aligned.toTranscriptome.out.bam",
        progress_log = "output/align_wasp/{sample}_Log.progress.out",
        tmp_dir = "output/align_wasp/{sample}_STARtmp",
        final_log = "output/align_wasp/{sample}_Log.final.out",
        sj_tab = "output/align_wasp/{sample}_SJ.out.tab"

    threads: 8
    
    log:
        err = 'output/logs/align_wasp_{sample}.err',
        out = 'output/logs/align_wasp_{sample}.out',
    
    params:
        index = config['star'],
        sjdb = config['starsjdb'],
        starVer=config['starVers'],
        dir_align = "output/align_wasp/{sample}_",
        vcf = lambda wildcards: generate_vcf_path(wildcards.sample)
    
    shell:
        """
        ml star/{params.starVer};
        mkdir -p output/align_wasp
        
        STAR --genomeDir {params.index} \
                --runThreadN {threads} \
                --sjdbFileChrStartEnd {params.sjdb} \
                --outFileNamePrefix {params.dir_align} \
                --quantMode TranscriptomeSAM \
                --outSAMstrandField intronMotif \
                --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --readFilesIn {input.R1} {input.R2} \
                --outFilterType BySJout \
                --outFilterMultimapNmax 20 \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --outFilterMismatchNmax 999 \
                --outFilterMismatchNoverReadLmax 0.04 \
                --alignIntronMin 20 \
                --alignIntronMax 1000000 \
                --alignMatesGapMax 1000000 \
                --waspOutputMode SAMtag \
                --varVCFfile {params.vcf} 1> {log.out} 2> {log.err}
        """

rule index:
    input:
        rules.align.output.bam
    output:
        bai = 'output/align_wasp/{sample}_Aligned.sortedByCoord.out.bam.bai'
    threads: 8
    params:
        samtoolsVersion = config['samtoolsVers']
    log:
        out = 'output/logs/index_wasp_{sample}.out',
        err = 'output/logs/index_wasp_{sample}.err'
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        samtools index -@ {threads} {input} {output} 1> {log.out} 2> {log.err}
        """


