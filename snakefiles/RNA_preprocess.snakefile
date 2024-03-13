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
        #[expand("output/fastq/{sampleName}_{read}.fastq.gz", sampleName=key,read=['R1','R2']) for key in read1],
        #[expand("output/QC/{sampleName}_{read}_fastqc.{ext}", sampleName=key, read=['R1', 'R2'], ext=['zip', 'html']) for key in read1],
        #[expand('output/trim/{sampleName}_R1{ext}', sampleName=key, ext=['.fastq.gz_trimming_report.txt', '_val_1.fq.gz']) for key in read1],
        #[expand('output/trim/{sampleName}_R2{ext}', sampleName=key, ext=['.fastq.gz_trimming_report.txt', '_val_2.fq.gz']) for key in read1],
        #[expand("output/{condition}_samples.txt", condition=['CTL','FNF','OA'])],
        #"output/geno/chQTL_samplesubset_103.vcf.gz",
        #[expand("output/geno/{condition}_geno/{condition}_matched.txt", condition=['pbs','fnf','oa'])],
        #[expand("output/geno/{condition}_geno/{condition}.{ext}", condition=['pbs','fnf','oa'],ext=["rename_wCHR.vcf.gz","rename.vcf.gz","rename_wCHR.vcf"])],
        #expand("output/align/{sampleName}_Aligned.sortedByCoord.out.bam", sampleName=read1.keys()),
        #expand("output/align/{sampleName}_Log.out", sampleName=read1.keys()),
        #expand("output/align/{sampleName}_Aligned.toTranscriptome.out.bam",sampleName=read1.keys()),
        #expand("output/align/{sampleName}_Log.progress.out", sampleName=read1.keys()),
        #expand("output/align/{sampleName}__STARtmp", sampleName=read1.keys()),
        #expand("output/align/{sampleName}_Log.final.out", sampleName=read1.keys()),
        #expand("output/align/{sampleName}_SJ.out.tab", sampleName=read1.keys()),
        #expand('output/align/{sampleName}_Aligned.sortedByCoord.out.bam.bai',sampleName=read1.keys()),        
        [expand('output/QC/{sampleName}_verifybamid.{ext}', sampleName=key,ext=['selfSM','selfRG','bestRG','bestSM','depthRG','depthSM'])for key in read1],
        [expand('output/align_wasp/{sampleName}.Aligned.sortedByCoord.WASP.{ext}', sampleName=key,ext=['bam','bam.bai'])for key in read1],
        [expand('output/quant/{sampleName}',sampleName=key) for key in read1]
        

include: '/work/users/s/e/seyoun/CQTL_sQTL/snakefiles/VCFpreprocess.snakefile'

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

rule reheader:
    input:
        VCF_in=rules.create_vcf.output.vcf,
        VCFOA_in=config['vcf_oa'],
        matched_pbs=rules.processVCF.output.pbs_matched,
        matched_fnf=rules.processVCF.output.fnf_matched,
        matched_oa=rules.processVCF.output.oa_matched
    output:
        pbs_vcf="output/geno/pbs_geno/pbs.rename_wCHR.vcf.gz",
        fnf_vcf="output/geno/fnf_geno/fnf.rename_wCHR.vcf.gz",
        oa_vcf="output/geno/oa_geno/oa.rename_wCHR.vcf.gz",
        pbs_temp="output/geno/pbs_geno/pbs.rename.vcf.gz",
        fnf_temp="output/geno/fnf_geno/fnf.rename.vcf.gz",
        oa_temp="output/geno/oa_geno/oa.rename.vcf.gz",
        nozip_pbs="output/geno/pbs_geno/pbs.rename_wCHR.vcf",
        nozip_fnf="output/geno/fnf_geno/fnf.rename_wCHR.vcf",
        nozip_oa="output/geno/oa_geno/oa.rename_wCHR.vcf"

    params:
        samtoolsVer=config['samtoolsVers'],
        plinkVer=config['plinkVers'],
        pbs_dir="output/geno/pbs_geno",
        fnf_dir="output/geno/fnf_geno",
        oa_dir="output/geno/oa_geno"
    log:
        pbsout="output/logs/PBS_rename.out",
        pbserr="output/logs/PBS_rename.err",
        fnfout="output/logs/FNF_rename.out",
        fnferr="output/logs/FNF_rename.err",
        oaout="output/logs/OA_rename.out",
        oaerr="output/logs/OA_rename.err"


    shell:
        """
        ml samtools/{params.samtoolsVer}
        ml plink/{params.plinkVer}

        sh /work/users/s/e/seyoun/CQTL_sQTL/scripts/VCF_rename.sh {input.VCF_in} {input.matched_pbs} {output.pbs_temp} {output.pbs_vcf} 1> {log.pbsout} 2> {log.pbserr}

        sh /work/users/s/e/seyoun/CQTL_sQTL/scripts/VCF_rename.sh {input.VCF_in} {input.matched_fnf} {output.fnf_temp} {output.fnf_vcf} 1> {log.fnfout} 2> {log.fnferr}

        sh /work/users/s/e/seyoun/CQTL_sQTL/scripts/VCF_rename.sh {input.VCFOA_in} {input.matched_oa} {output.oa_temp} {output.oa_vcf} 1> {log.oaout} 2> {log.oaerr}

        gunzip -c {output.pbs_vcf} > {output.nozip_pbs}
        gunzip -c {output.fnf_vcf} > {output.nozip_fnf}
        gunzip -c {output.oa_vcf} > {output.nozip_oa}


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
        R2 = rules.trim.output.trim2,
        vcf_default= rules.reheader.output.nozip_pbs

    output:
        bam = "output/align/{sampleName}_Aligned.sortedByCoord.out.bam",
        log_out = "output/align/{sampleName}_Log.out",
        #transcriptome_bam = "output/align/{sampleName}_Aligned.toTranscriptome.out.bam",
        progress_log = "output/align/{sampleName}_Log.progress.out",
        #tmp_dir = "output/align/{sampleName}__STARtmp",
        final_log = "output/align/{sampleName}_Log.final.out",
        #sj_tab = "output/align/{sampleName}_SJ.out.tab"

    threads: 8

    log:
        err = 'output/logs/align_{sampleName}.err',
        out = 'output/logs/align_{sampleName}.out'
    params:
        index = config['star'],
        sjdb = config['starsjdb'],
        starVer=config['starVers'],
        dir_align = "output/align/{sampleName}_"
    shell:
        """
        ml star/{params.starVer};
        mkdir -p output/align
        echo '{input.vcf_default}'

        sampleName=$(basename {input.R1} | cut -d '_' -f 2)
        echo "$sampleName"

        if [[ $sampleName == *"CTL"* ]]; then
            vcf_path="output/geno/pbs_geno/pbs.rename_wCHR.vcf"
        elif [[ $sampleName == *"FNF"* ]]; then
            vcf_path="output/geno/fnf_geno/fnf.rename_wCHR.vcf"
        elif [[ $sampleName == *"OA"* ]]; then
            vcf_path="output/geno/oa_geno/oa.rename_wCHR.vcf"
        else
        echo "Error: Unrecognized sampleName $sampleName"
        exit 1
        fi
        
        echo "$vcf_path"

        STAR --genomeDir {params.index} \
                --runThreadN {threads} \
                --sjdbFileChrStartEnd {params.sjdb} \
                --outFileNamePrefix {params.dir_align} \
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
                --varVCFfile $vcf_path 1> {log.out} 2> {log.err}
        """

rule index:
    input:
        rules.align.output.bam
    output:
        temp='output/align/{sampleName}_Aligned.sortedByCoord.out.bam.bai',
        bam= 'output/align/{sampleName}_sorted.bam',
        bai = 'output/align/{sampleName}_sorted.bai'
    threads: 8
    params:
        samtoolsVersion = config['samtoolsVers']
    log:
        out = 'output/logs/index_{sampleName}.out',
        err = 'output/logs/index_{sampleName}.err'
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        samtools index -@ {threads} {input} {output.temp} 
        samtools view -b {input} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY |samtools sort -o {output.bam} 
        samtools index -@ {threads} {output.bam} {output.bai} 2> {log.err}
        """

rule addReadGroups:
    input:
        bam = rules.index.output.bam,
        bai = rules.index.output.bai
    output:
        bam = 'output/align/{sampleName}_sorted.RG.bam',
        bai = 'output/align/{sampleName}_sorted.RG.bai'
    params:
        picardVer =  config['picardVers'],
        samtoolsVersion = config['samtoolsVers']
    log:
        err1 = 'output/logs/addReadGroups_{sampleName}.err',
        err2 = 'output/logs/addReadGroups_index_{sampleName}.err'
    shell:
        """
        module load picard/{params.picardVer}
        module load samtools/{params.samtoolsVersion}

        # Assuming the input string
        trimmed_string={input.bam}
        # Remove the prefix
        
        trimmed_string=${{trimmed_string#output/align/}}
        # Extracting sampleName and phase
        IFS="_" read -r -a array <<< "${{trimmed_string}}"
        donor="${{array[0]}}"
        
        picard AddOrReplaceReadGroups -I {input.bam} -O {output.bam} --RGSM $donor --RGPL ILLUMINA --RGLB lib1 --RGPU unit1 2> {log.err1}
        # Index
        samtools index -@ {threads} {output.bam} {output.bai} 2> {log.err2}
        """

rule verifybamid:
    input:
        bam = rules.addReadGroups.output.bam,
        bai = rules.addReadGroups.output.bai
    output:
        selfSM = 'output/QC/{sampleName}_verifybamid.selfSM',
        selfRG = 'output/QC/{sampleName}_verifybamid.selfRG',
        bestRG = 'output/QC/{sampleName}_verifybamid.bestRG',
        bestSM = 'output/QC/{sampleName}_verifybamid.bestSM',
        depthRG = 'output/QC/{sampleName}_verifybamid.depthRG',
        depthSM = 'output/QC/{sampleName}_verifybamid.depthSM'
    params:
        verifybamid = config['verifybamid'],
        out_dir= 'output/QC/{sampleName}_verifybamid'
    benchmark:
        'output/benchmarks/{sampleName}__verifybamid.tsv'
    log:
        err = 'output/logs/verifybamid_{sampleName}.err'
    shell:
        """
        sample=$(basename {input.bam} | cut -d '_' -f 2)
        echo "$sample"

        if [[ $sample == *"CTL"* ]]; then
            vcf_path="output/geno/pbs_geno/pbs.rename_wCHR.vcf"
        elif [[ $sample == *"FNF"* ]]; then
            vcf_path="output/geno/fnf_geno/fnf.rename_wCHR.vcf"
        elif [[ $sample == *"OA"* ]]; then
            vcf_path="output/geno/oa_geno/oa.rename_wCHR.vcf"
        else
            echo "Error: Unrecognized sampleName $sampleName"
            exit 1
        fi

        echo "$vcf_path"
        
        {params.verifybamid} --vcf $vcf_path --bam {input.bam} --bai {input.bai} --best --out {params.out_dir} 2> {log.err}
        """

#Filter tagged WASP reads
rule WASPfilter:
    input:
        R = rules.index.output.bam
    output:
        bam = 'output/align_wasp/{sampleName}.Aligned.sortedByCoord.WASP.bam',
        bai = 'output/align_wasp/{sampleName}.Aligned.sortedByCoord.WASP.bam.bai'
    threads: 2
    log:
        err1 = 'output/logs/WASPfilter_{sampleName}_grep.err',
        err2 = 'output/logs/WASPfilter_{sampleName}_samtoolsView.err',
        err3 = 'output/logs/WASPfilter_{sampleName}_samtoolsIndex.err'
    params:
        samtoolsVersion = config['samtoolsVers']
    shell:
        """
        module load samtools/{params.samtoolsVersion}
        mkdir -p output/align_wasp
        # Add header
        samtools view -H {input.R} > output/align_wasp/{wildcards.sampleName}.Aligned.sortedByCoord.WASP.sam
        
        # Grep for WASP-passing reads
        samtools view {input.R} | grep 'vW:i:1' >> output/align_wasp/{wildcards.sampleName}.Aligned.sortedByCoord.WASP.sam 2> {log.err1}

        # Compress and index
        samtools view -bS output/align_wasp/{wildcards.sampleName}.Aligned.sortedByCoord.WASP.sam > {output.bam} 2> {log.err2}
        samtools index -@ {threads} {output.bam} {output.bai} 2> {log.err3}
        """


rule quant:
    input:
        trim1 = rules.trim.output.trim1,
        trim2 = rules.trim.output.trim2
    output:
        'output/quant/{sampleName}'
    params:
        salmonVer = config['salmonVers'],
        index=config['salmon'],
        gcFlag = config['gcBias'],
        seqFlag = config['seqBias']
    log:
        out='output/logs/salmon_{sampleName}.out',
        err='output/logs/salmon_{sampleName}.err'
    shell:
        """
        ml salmon/{params.salmonVer}
        mkdir -p output/quant
        salmon quant -i {params.index} -l A -1 {input.trim1} -2 {input.trim2} -o {output} --seqBias --gcBias --validateMappings
        """

