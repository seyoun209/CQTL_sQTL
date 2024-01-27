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
rule all:
    input:
        #[expand("output/fastq/{sampleName}_{read}.fastq.gz", sampleName=key,read=['R1','R2']) for key in read1]
        [expand("output/QC/{sampleName}_{read}_fastqc.{ext}", sampleName=key, read=['R1', 'R2'], ext=['zip', 'html']) for key in read1],
        [expand('output/trim/{sampleName}_{R}{ext}',sampleName=key, R=['R1', 'R2'],ext=['_trimmed.fq.gz', '.fastq.gz_trimming_report.txt']) for key in read1]
        #[expand("output/fastq_unzip/{sampleName}_{read}.fastq", sampleName=key,read=['R1','R2']) for key in read1]
        # [expand("output/align/{sampleName}_Aligned.sortedByCoord.out.bam", sampleName=key) for key in read1]
        #[expand("output/align/{sampleName}{ext}", sampleName=key,ext=align_output) for key in read1]
        #[expand("output/align/{case}.txt",case=['control','case'])]
        #"output/rmats"
        #("output/QC/{name}_multiqc_report.html").format(name=runName)
        #expand("output/align/{study}.txt",study=['control','case'])

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
        trim1 = 'output/trim/{sampleName}_R1_trimmed.fq.gz',
        trim2 = 'output/trim/{sampleName}_R2_trimmed.fq.gz',
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
        module load python/3.6.6
        module load pigz
        mkdir -p output/trim
        trim_galore -o output/trim --cores {threads} --path_to_cutadapt /nas/longleaf/apps/cutadapt/2.9/venv/bin/cutadapt --paired {input.R1} {input.R2} 2> {log.err}
        """
