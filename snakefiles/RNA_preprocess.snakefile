#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os
from utils.namer import namer
#from snakemake.io import *
#from workflows.utils.namer import namer

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep='\t')
#samples = pd.read_csv("samplesheet_NA_removed.txt",sep='\t')
## Convert all columns to strings
samples = samples.astype(str)

## Concatenate the sequencing directory to Read1 and Read2 for full paths
samples['Read1'] = samples[['Sequencing_Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)
samples['Read2'] = samples[['Sequencing_Directory', 'Read2']].apply(lambda row: os.path.join(*row), axis=1)

## Concatenate columns to identify which groups to run (i.e. Seq_Rep will be run together)
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)
#mergeBy = ['Proj','Donor','Condition','Time']
#samples['mn'] = samples[mergeBy].agg('_'.join, axis=1)

## Group by mn and extract Read1 & Read2
read1 = samples.groupby('mn')['Read1'].apply(list).to_dict()
read2 = samples.groupby('mn')['Read2'].apply(list).to_dict()

## Set run summary name using helper script
runName = namer(samples, config['mergeBy'])
#runName = namer(samples, samples[mergeBy])

#align output format
align_output = ["_STARgenome",".Aligned.sortedByCoord.out.bam","Aligned.toTranscriptome.out.bam","Log.final.out","Log.out","Log.progress.out","SJ.out.tab"]
#[".Aligned.sortedByCoord.out.bam","Aligned.toTranscriptome.out.bam"]
#
## Group by mn to update samplesheet 
# reads = [expand("{sampleName}_{read}.fastq.gz", sampleName=key,read=['R1','R2']) for key in read1]
# namePD = pd.DataFrame.from_dict(dict(zip(read1,reads)), orient='index')
# namePD.index.name="mn"
# newsamplesheet = samples.merge(namePD,how='left',on='mn')
# newsamplesheet['new_directory'] = os.getcwd()+'/output/fastq'
# newsamplesheet.rename(columns={0:'New_read1',1:'New_read2'},inplace=True)
# newsamplesheet.to_csv("output/samplesheetupdated.txt", sep="\t",index=False)

## Define actions on success
onsuccess:
	## Success message
	print("RNA_Splicing completed successfully! Wahoo!")

##### Define rules #####
include: '/work/users/s/e/seyoun/CQTL_sQTL/snakefiles/VCFpreprocess.snakefile'
rule all:
    input:
        #[expand("output/fastq/{sampleName}_{read}.fastq.gz", sampleName=key,read=['R1','R2']) for key in read1],
        #[expand("output/QC/{sampleName}_{read}_fastqc.{ext}", sampleName=key, read=['R1', 'R2'], ext=['zip', 'html']) for key in read1],
        #[expand('output/trim/{sampleName}_R1{ext}', sampleName=key, ext=['.fastq.gz_trimming_report.txt', '_val_1.fq.gz']) for key in read1],
        #[expand('output/trim/{sampleName}_R2{ext}', sampleName=key, ext=['.fastq.gz_trimming_report.txt', '_val_2.fq.gz']) for key in read1],
        [expand("output/{condition}_samples.txt", condition=['CTL','FNF','OA'])],
        "output/geno/chQTL_samplesubset_103.vcf.gz",
        [expand("output/geno/{condition}_geno/{condition}_matched.txt", condition=['pbs','fnf','oa'])],
        [expand("output/geno/{condition}_geno/{condition}.{ext}", condition=['pbs','fnf','oa'],ext=["rename_wCHR.vcf.gz","rename.vcf.gz"])],
        #expand("output/align/{sampleName}_Aligned.sortedByCoord.out.bam", sampleName=read1.keys()),
        #expand("output/align/{sampleName}_Log.out", sampleName=read1.keys()),
        #expand("output/align/{sampleName}_Aligned.toTranscriptome.out.bam",sampleName=read1.keys()),
        #expand("output/align/{sampleName}_Log.progress.out", sampleName=read1.keys()),
        #expand("output/align/{sampleName}__STARtmp", sampleName=read1.keys()),
        #expand("output/align/{sampleName}_Log.final.out", sampleName=read1.keys()),
        #expand("output/align/{sampleName}_SJ.out.tab", sampleName=read1.keys()),
        expand('output/align/{sampleName}_Aligned.sortedByCoord.out.bam.bai',sampleName=read1.keys())




rule catReads:
    input:
        R1 = lambda wildcards: read1.get(wildcards.sampleName),
        R2 = lambda wildcards: read2.get(wildcards.sampleName)
    output:
        R1 = 'output/fastq/{sampleName}_R1.fastq.gz',
        R2 = 'output/fastq/{sampleName}_R2.fastq.gz'
    benchmark:
        'output/benchmarks/{sampleName}_catReads.tsv'
    log:
        errR1 = 'output/logs/{sampleName}_R1_catReads.err',
        errR2 = 'output/logs/{sampleName}_R2_catReads.err'
    shell:
        """
        cat {input.R1} > {output.R1} 2> {log.errR1}
        cat {input.R2} > {output.R2} 2> {log.errR2}
        """

rule fastqc:
	input:
		#R1 = lambda wildcards: read1.get(wildcards.sampleName),
		#R2 = lambda wildcards: read2.get(wildcards.sampleName)
		QC1 = rules.catReads.output.R1,
		QC2 = rules.catReads.output.R2
	output:
		zip1 = "output/QC/{sampleName}_R1_fastqc.zip", # temp
		zip2 = "output/QC/{sampleName}_R2_fastqc.zip", # temp
		html1 = "output/QC/{sampleName}_R1_fastqc.html", # temp
		html2 = "output/QC/{sampleName}_R2_fastqc.html" # temp
	log:
		err = 'output/logs/fastqc_{sampleName}.err',
		out = 'output/logs/fastqc_{sampleName}.out'
	params:
		dir = "output/QC",
		version = config['fastqcVers']
	benchmark: 
		'output/benchmarks/fastqc_{sampleName}.tsv'
	shell:
		"""
		module load fastqc/{params.version};
		fastqc -o {params.dir} {input.QC1} {input.QC2} 1> {log.out} 2> {log.err};
		"""          
rule trim:
    input:
        R1 = rules.catReads.output.R1,
        R2 = rules.catReads.output.R2
    output:
        trim1 = 'output/trim/{sampleName}_R1_val_1.fq.gz',
        trim2 = 'output/trim/{sampleName}_R2_val_2.fq.gz',
        report1 = 'output/trim/{sampleName}_R1.fastq.gz_trimming_report.txt',
        report2 = 'output/trim/{sampleName}_R2.fastq.gz_trimming_report.txt'
    threads: 4
    params:
        version = config['trim_galore']
    log:
        err = 'output/logs/trim_{sampleName}.err',
        out = 'output/logs/trim_{sampleName}.out'
    shell:
        """
        module load trim_galore/{params.version}
        module load python/3.9.6
        module load pigz
        mkdir -p output/trim
        trim_galore -o output/trim --cores {threads} --paired {input.R1} {input.R2} 2> {log.err}
        """

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
        R1 = rules.trim.output.trim1,
        R2 = rules.trim.output.trim2
    output:
        bam = "output/align/{sampleName}_Aligned.sortedByCoord.out.bam",
        log_out = "output/align/{sampleName}_Log.out",
        transcriptome_bam = "output/align/{sampleName}_Aligned.toTranscriptome.out.bam",
        progress_log = "output/align/{sampleName}_Log.progress.out",
        tmp_dir = "output/align/{sampleName}__STARtmp",
        final_log = "output/align/{sampleName}_Log.final.out",
        sj_tab = "output/align/{sampleName}_SJ.out.tab"

    threads: 8

    log:
        err = 'output/logs/align_{sampleName}.err',
        out = 'output/logs/align_{sampleName}.out',
    params:
        index = config['star'],
        sjdb = config['starsjdb'],
        starVer=config['starVers'],
        dir_align = "output/align/{sampleName}_",
        vcf = generate_vcf_path(rules.trim.output.trim1)
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
        bai = 'output/align/{sampleName}_Aligned.sortedByCoord.out.bam.bai'
    threads: 8
    params:
        samtoolsVersion = config['samtoolsVers']
    log:
        out = 'output/logs/index_{sampleName}.out',
        err = 'output/logs/index_{sampleName}.err'
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        samtools index -@ {threads} {input} {output} 1> {log.out} 2> {log.err}
        """


