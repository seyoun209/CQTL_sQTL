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



## Define actions on success
onsuccess:
	## Success message
	print("RNApipe completed successfully! Wahoo!")

##### Define rules #####        
rule all:
    input:
        expand('output/rmats_{case}',case=['fnf','oa']),
        expand('output/junc/{sampleName}.junc',sampleName=bamfile.keys()),
        expand('output/junc_wasp/{sampleName}_wasp.junc',sampleName=bamfile.keys()),
        #'output/clu_fnf/ctlvsfnf_perind_numers.counts.gz',
        #'output/clu_oa/ctlvsoa_perind_numers.counts.gz',
        #'output/clu_fnf_wasp/ctlvsfnf_perind_numers.counts.gz',
        #'output/clu_oa_wasp/ctlvsoa_perind_numers.counts.gz'
        

rule samplelist:
    input:
        expand("output/align/{sampleName}_sorted.bam", sampleName=key) for key in bamfile
    output:
        txt1= "output/align/ctl.txt",
        txt2= "output/align/fnf.txt",
        txt3= "output/align/oa.txt"
    log:
        err1 = 'output/logs/spList_ctl.err',
	err2 = 'output/logs/spList_fnf.err',
        err3 = 'output/logs/spList_oa.err'
    shell:
        """
        find ./output/align/*CTL*  -name *_sorted.bam* |paste -sd, -  > {output.txt1} 2> {log.err1}
        find ./output/align/*FNF*  -name *_sorted.bam* |paste -sd, -  > {output.txt2} 2> {log.err2}
        find ./output/align/*OA*  -name *_sorted.bam* |paste -sd, -  > {output.txt3} 2> {log.err3}
        """
rule rmates:
    input:
        ctrList = rules.samplelist.output.txt1,
        fnfList = rules.samplelist.output.txt2,
        oaList = rules.samplelist.output.txt3
    output:
        dir_fnf= 'output/rmats_fnf',
        dir_oa= 'output/rmats_oa'
    log:
        err1 = 'output/rmats_fnf/rmats_ctl_fnf',
        err2 = 'output/rmats_oa/rmats_ctl_oa'
    params:
        rmatsVersion = config['rmats_turbo'],
        gtf = config['gtf'],
        strand = config['t'],
        rlength= config['readlength'],
    threads: 4
    shell:
        """
        module load rmats-turbo/{params.rmatsVersion}

        run_rmats --b1 {input.ctrList} --b2 {input.fnfList} --gtf {params.gtf}  -t {params.strand} --readLength {params.rlength} --nthread {threads} --od {output.dir_fnf}  --tmp {log.err1}

        run_rmats --b1 {input.ctrList} --b2 {input.oaList} --gtf {params.gtf} -t {params.strand} --readLength {params.rlength} --nthread {threads} --od {output.dir_oa} --tmp {log.err2}
        """


rule junc:
    input:
        bam1= lambda wildcards: bamfile.get(wildcards.sampleName),
        bam2= lambda wildcards: bamfilewasp.get(wildcards.sampleName)
    output:
        junc = 'output/junc/{sampleName}.junc',
        junc_wasp = 'output/junc_wasp/{sampleName}_wasp.junc'
    params:
        regtoolVer = config['regtools']
    log:
        err = 'output/logs/regtools_{sampleName}.err',
        err_wasp = 'output/logs/regtools_wasp_{sampleName}.err'
    shell:
        """
        mkdir -p output/junc
        mkdir -p output/junc_wasp
        
        {params.regtoolVer} junctions extract -a 8 -m 50 -M 500000 -s XS {input.bam1} -o {output.junc} 2> {log.err}
        {params.regtoolVer} junctions extract -a 8 -m 50 -M 500000 -s XS {input.bam2} -o {output.junc_wasp}  2> {log.err_wasp}
        """











