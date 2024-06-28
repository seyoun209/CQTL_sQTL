#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os
from utils.namer import namer
#from snakemake.io import *
#from snakefiles.utils.namer import namer


# Read in samplesheet
samples = pd.read_csv(config["update_samplesheet"], sep='\t')
samples = samples.astype(str)

# Omit samples
samples = samples[~samples['Donor'].str.contains('|'.join(config['samples_to_omit']))]

# Create mn column
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)


# Group by mn and bam
bamfile = samples.groupby('mn')['alignbam'].apply(list).to_dict()
bamfilewasp= samples.groupby('mn')['alignbam_wasp'].apply(list).to_dict()

## Build dictionary of merged BAM files
mergeSample = samples.groupby('Condition')['mn'].apply(list).to_dict()

# Unique sample names in dictionary
#mergeSamples_dedup = dict()
#for key in mergeSample:
#    mergeSamples_dedup[key] = list(set(mergeSample[key]))

##### Define rules #####
rule all:
    input:
        #outfiles
        [expand("output/signals/01.signal_fnf/{sampleName}.bw", sampleName=key) for key in bamfile],
        [expand("output/signals/02.signal_fnf_wasp/{sampleName}.bw", sampleName=key) for key in bamfilewasp],
        [expand("output/signals/03.mergeAlign/{cond}_{ext}", cond=key, ext=['sorted.bam', 'sorted.bam.bai', 'stats.txt']) for key in mergeSample],
        [expand("output/signals/04.mergeAlign_wasp/{cond}_{ext}", cond=key, ext=['sorted.bam', 'sorted.bam.bai', 'stats.txt']) for key in mergeSample],
        [expand("output/signals/05.mergesignal/{cond}.bw", cond=key) for key in mergeSample],
        expand("output/signals/strand/{cond}_fwd.bw",cond=mergeSample.keys()),
        expand("output/signals/strand/{cond}_rev.bw",cond=mergeSample.keys()),
        expand("output/signals/merged_norm/{cond}_norm.bw",cond=mergeSample.keys())
        

rule signal:
    input:
        bam= lambda wildcards: bamfile.get(wildcards.sampleName),
        bam_wasp= lambda wildcards: bamfilewasp.get(wildcards.sampleName)
    output:
        signal = "output/signals/01.signal_fnf/{sampleName}.bw",
        signal_wasp = "output/signals/02.signal_fnf_wasp/{sampleName}.bw"
    log:
        err1='output/logs/signal_{sampleName}.err',
        err2='output/logs/signal_wasp_{sampleName}.err'

    params:
        deeptools_ver = config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p output/signals/01.signal_fnf
        mkdir -p output/signals/02.signal_fnf_wasp

        bamCoverage -b {input.bam} -o {output.signal} > {log.err1} 2>&1
        bamCoverage -b {input.bam_wasp} -o {output.signal_wasp} > {log.err2} 2>&1

        """

rule mergeAlign:
    input:
        bam_files = lambda wildcards: ["output/align/{sampleName}_sorted.bam".format(sampleName=value) for value in mergeSample[wildcards.cond]],
        bam_wasp_files = lambda wildcards: ["output/align_wasp/{sampleName}.Aligned.sortedByCoord.WASP.bam".format(sampleName=value) for value in mergeSample[wildcards.cond]]
    
    output:
        bam = "output/signals/03.mergeAlign/{cond}_sorted.bam",
	bai = "output/signals/03.mergeAlign/{cond}_sorted.bam.bai",
	stats = "output/signals/03.mergeAlign/{cond}_stats.txt",
        bam_wasp = "output/signals/04.mergeAlign_wasp/{cond}_sorted.bam",
        bai_wasp = "output/signals/04.mergeAlign_wasp/{cond}_sorted.bam.bai",
        stats_wasp = "output/signals/04.mergeAlign_wasp/{cond}_stats.txt"
    log:
        err1 = 'output/logs/mergeAlign_{cond}.err',
        err2 = 'output/logs/mergeAlign_wasp_{cond}.err'
    params:
        samtoolsVersion = config['samtoolsVers']
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        mkdir -p output/signals/03.mergeAlign
        mkdir  -p output/signals/04.mergeAlign_wasp

        samtools merge {output.bam} {input.bam_files} >> {log.err1} 2>&1
	samtools flagstat {output.bam} > {output.stats} >> {log.err1} 2>&1
	samtools index {output.bam} >> {log.err1} 2>&1

        samtools merge {output.bam_wasp} {input.bam_wasp_files} >> {log.err2} 2>&1
        samtools flagstat {output.bam_wasp} > {output.stats_wasp} >> {log.err2} 2>&1 
        samtools index {output.bam_wasp} >> {log.err2} 2>&1

        """

rule mergeSignal:
    input:
        bam = rules.mergeAlign.output.bam
    output:
        signal = "output/signals/05.mergesignal/{cond}.bw"
    log:
        err1 = 'output/logs/mergeSignal_{cond}.err',
    params:
        deeptools_ver = config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p output/signals/05.mergesignal
        
        bamCoverage -b {input.bam} -o {output.signal} > {log.err1} 2>&1
        
        
        """

rule mergeForwadSignal:
    input:
        bam = rules.mergeAlign.output.bam
    output:
        signal_fw = "output/signals/strand/{cond}_fwd.bw"
    log:
        err1 = 'output/logs/mergeSignal_forward_{cond}.err'
    params:
        deeptools_ver = config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p output/signals/strand
        
            
        bamCoverage --filterRNAstrand forward bamCoverage --normalizeUsing BPM --binSize 10 --effectiveGenomeSize 2862010578 -b {input.bam} -o {output.signal_fw} > {log.err1} 2>&1
        

        """

rule mergeReverseSignal:
    input:
        bam = rules.mergeAlign.output.bam
    output:
        signal_rv = "output/signals/strand/{cond}_rev.bw"
    log:
        err1 = 'output/logs/mergeSignal_reverse_{cond}.err'
    params:
        deeptools_ver = config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p output/signals/strand
        

        bamCoverage --filterRNAstrand reverse bamCoverage --normalizeUsing BPM --binSize 10 --effectiveGenomeSize 2862010578 -b {input.bam} -o {output.signal_rv} > {log.err1} 2>&1
        
            
        """

rule mergeSignal_norm:
    input:
        bam= rules.mergeAlign.output.bam
    output:
        signal = "output/signals/merged_norm/{cond}_norm.bw"
    log:
        err1 = 'output/logs/mergeSignal_normalized_{cond}.err'
    params:
        deeptools_ver = config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p output/signals/merged_norm

        bamCoverage --normalizeUsing BPM --binSize 10 --effectiveGenomeSize 2862010578 --bam {input.bam} -o {output.signal} > {log.err1} 2>&1
        """

