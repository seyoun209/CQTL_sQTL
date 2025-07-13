#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os
from utils.namer import namer

# Read in samplesheet and filter
samples = pd.read_csv(config["update_samplesheet"], sep='\t').astype(str)
samples = samples[~samples['Donor'].str.contains('|'.join(config['samples_to_omit']))]
samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)
bamfile = samples.groupby('mn')['alignbam'].apply(list).to_dict()

onsuccess:
    print("RNApipe completed successfully! Wahoo!")

##### Define rules #####        
rule all:
    input:
        expand('output/rmats_{case}', case=['fnf','oa']),
        expand('output/junc/{sampleName}.junc', sampleName=bamfile.keys())
        # Add further outputs for downstream steps as needed

rule samplelist:
    input:
        expand("output/align/{sampleName}_sorted.bam", sampleName=key) for key in bamfile
    output:
        txt1="output/align/ctl.txt",
        txt2="output/align/fnf.txt",
        txt3="output/align/oa.txt"
    log:
        err1='output/logs/spList_ctl.err',
        err2='output/logs/spList_fnf.err',
        err3='output/logs/spList_oa.err'
    shell:
        """
        find ./output/align/*CTL*  -name *_sorted.bam* | paste -sd, -  > {output.txt1} 2> {log.err1}
        find ./output/align/*FNF*  -name *_sorted.bam* | paste -sd, -  > {output.txt2} 2> {log.err2}
        find ./output/align/*OA*   -name *_sorted.bam* | paste -sd, -  > {output.txt3} 2> {log.err3}
        """

rule extract_junctions:
    input:
        bam="output/align/{sampleName}_sorted.bam"
    output:
        junc="output/junc/{sampleName}.junc"
    log:
        out="output/logs/extract_junctions_{sampleName}.out",
        err="output/logs/extract_junctions_{sampleName}.err"
    params:
        regtools=config['regtools']
    shell:
        """
        mkdir -p output/junc
        {params.regtools} junctions extract -a 8 -m 50 -M 500000 {input.bam} -o {output.junc} > {log.out} 2> {log.err}
        """

rule rmates:
    input:
        ctrList=rules.samplelist.output.txt1,
        fnfList=rules.samplelist.output.txt2,
        oaList=rules.samplelist.output.txt3
    output:
        dir_fnf='output/rmats_fnf',
        dir_oa='output/rmats_oa'
    log:
        err1='output/rmats_fnf/rmats_ctl_fnf',
        err2='output/rmats_oa/rmats_ctl_oa'
    params:
        rmatsVersion=config['rmats_turbo'],
        gtf=config['gtf'],
        strand=config['t'],
        rlength=config['readlength']
    threads: 4
    shell:
        """
        mkdir -p output/rmats_fnf output/rmats_oa
        # Example: run rMATS for FNF vs CTL and OA vs CTL
        python {params.rmatsVersion} --b1 {input.ctrList} --b2 {input.fnfList} \
            --gtf {params.gtf} --readLength {params.rlength} --variable-read-length \
            --nthread {threads} --od output/rmats_fnf --tmp output/rmats_fnf/tmp --libType fr-unstranded > {log.err1} 2>&1

        python {params.rmatsVersion} --b1 {input.ctrList} --b2 {input.oaList} \
            --gtf {params.gtf} --readLength {params.rlength} --variable-read-length \
            --nthread {threads} --od output/rmats_oa --tmp output/rmats_oa/tmp --libType fr-unstranded > {log.err2} 2>&1
        """
rule junc:
    input:
        bam = "output/align/{sampleName}_sorted.bam"
    output:
        junc = "output/junc/{sampleName}.junc"
    log:
        out = "output/logs/junc_{sampleName}.out",
        err = "output/logs/junc_{sampleName}.err"
    params:
        regtools = config['regtools']
    shell:
        """
        mkdir -p output/junc
        {params.regtools} junctions extract -a 8 -m 50 -M 500000 {input.bam} -o {output.junc} > {log.out} 2> {log.err}
        """
