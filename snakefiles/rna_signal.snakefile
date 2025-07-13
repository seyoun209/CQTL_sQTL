#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os
from utils.namer import namer

# Read in samplesheet and convert all columns to strings
samples = pd.read_csv(config["update_samplesheet"], sep='\t').astype(str)

# Omit samples
samples = samples[~samples['Donor'].str.contains('|'.join(config['samples_to_omit']))]

# Create mn column
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)

# Group by mn and BAM
bamfile = samples.groupby('mn')['alignbam'].apply(list).to_dict()
mergeSample = samples.groupby('Condition')['mn'].apply(list).to_dict()

##### Define rules #####
rule all:
    input:
        [expand("output/signals/01.signal_fnf/{sampleName}.bw", sampleName=key) for key in bamfile],
        [expand("output/signals/03.mergeAlign/{cond}_{ext}", cond=key, ext=['sorted.bam', 'sorted.bam.bai', 'stats.txt']) for key in mergeSample],
        [expand("output/signals/05.mergesignal/{cond}.bw", cond=key) for key in mergeSample],
        expand("output/signals/strand/{cond}_fwd.bw", cond=mergeSample.keys()),
        expand("output/signals/strand/{cond}_rev.bw", cond=mergeSample.keys()),
        expand("output/signals/merged_norm/{cond}_norm.bw", cond=mergeSample.keys())

rule signal:
    input:
        bam=lambda wildcards: bamfile.get(wildcards.sampleName)
    output:
        signal="output/signals/01.signal_fnf/{sampleName}.bw"
    log:
        err='output/logs/signal_{sampleName}.err'
    params:
        deeptools_ver=config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p output/signals/01.signal_fnf
        bamCoverage -b {input.bam} -o {output.signal} > {log.err} 2>&1
        """

rule mergeAlign:
    input:
        bam_files=lambda wildcards: [f"output/align/{mn}_Aligned.sortedByCoord.out.bam" for mn in mergeSample[wildcards.cond]]
    output:
        bam="output/signals/03.mergeAlign/{cond}_sorted.bam",
        bai="output/signals/03.mergeAlign/{cond}_sorted.bam.bai",
        stats="output/signals/03.mergeAlign/{cond}_stats.txt"
    log:
        err1='output/logs/mergeAlign_{cond}.err'
    params:
        samtools_ver=config['samtoolsVers']
    shell:
        """
        module load samtools/{params.samtools_ver}
        mkdir -p output/signals/03.mergeAlign
        samtools merge {output.bam} {input.bam_files} >> {log.err1} 2>&1
        samtools flagstat {output.bam} > {output.stats} >> {log.err1} 2>&1
        samtools index {output.bam} >> {log.err1} 2>&1
        """

rule mergeSignal:
    input:
        bam=rules.mergeAlign.output.bam
    output:
        signal="output/signals/05.mergesignal/{cond}.bw"
    log:
        err1='output/logs/mergeSignal_{cond}.err'
    params:
        deeptools_ver=config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p output/signals/05.mergesignal
        bamCoverage -b {input.bam} -o {output.signal} > {log.err1} 2>&1
        """

rule mergeForwadSignal:
    input:
        bam=rules.mergeAlign.output.bam
    output:
        signal_fw="output/signals/strand/{cond}_fwd.bw"
    log:
        err1='output/logs/mergeSignal_forward_{cond}.err'
    params:
        deeptools_ver=config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p output/signals/strand
        bamCoverage --filterRNAstrand forward --normalizeUsing BPM --binSize 10 --effectiveGenomeSize 2862010578 -b {input.bam} -o {output.signal_fw} > {log.err1} 2>&1
        """

rule mergeReverseSignal:
    input:
        bam=rules.mergeAlign.output.bam
    output:
        signal_rv="output/signals/strand/{cond}_rev.bw"
    log:
        err1='output/logs/mergeSignal_reverse_{cond}.err'
    params:
        deeptools_ver=config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p output/signals/strand
        bamCoverage --filterRNAstrand reverse --normalizeUsing BPM --binSize 10 --effectiveGenomeSize 2862010578 -b {input.bam} -o {output.signal_rv} > {log.err1} 2>&1
        """

rule mergeSignal_norm:
    input:
        bam=rules.mergeAlign.output.bam
    output:
        signal="output/signals/merged_norm/{cond}_norm.bw"
    log:
        err1='output/logs/mergeSignal_normalized_{cond}.err'
    params:
        deeptools_ver=config['deeptools']
    shell:
        """
        module load deeptools/{params.deeptools_ver}
        mkdir -p output/signals/merged_norm
        bamCoverage --normalizeUsing BPM --binSize 10 --effectiveGenomeSize 2862010578 -b {input.bam} -o {output.signal} > {log.err1} 2>&1
        """
