#!/usr/bin/env python3

CHR_VALUES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10","11", "12", "13", "14", "15", "16", "17", "18", "19",
              "20", "21", "22"]

rule all:
    input:
        "output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz",
        "output/geno/pbs_geno/pbs_rename_101_wCHR.vcf.gz",
        "output/geno/pbs_geno/pbs_rename_101.vcf.gz",
        "output/geno/fnf_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz",
        "output/geno/fnf_geno/fnf_rename_101_wCHR.vcf.gz",
        "output/geno/fnf_geno/fnf_rename_101.vcf.gz",
        "output/geno/oa_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz",
        "output/geno/oa_geno/oa_rename_101_wCHR.vcf.gz",
        "output/geno/oa_geno/oa_rename_101.vcf.gz",
        #"output/gtex_cluster/qtltools_prep/covariates_PC1",
        #"output/gtex_cluster_wasp/qtltools_prep/covariates_PC1"
        expand("output/gtex_cluster/qtltools_prep/covariates_PC{pc}",pc=range(1,21)),
        expand("output/gtex_cluster_wasp/qtltools_prep/covariates_PC{pc}",pc=range(1,21)),
        expand("output/01.qtltools_re/nominal_pbs/pc{pc}/chr{chr}.pbs.cis", pc=range(1, 21), chr=CHR_VALUES),
        expand("output/01.qtltools_re/nominal_fnf/pc{pc}/chr{chr}.fnf.cis", pc=range(1, 21), chr=CHR_VALUES),
        expand("output/01.qtltools_re/perm_pbs/pc{pc}/chr{chr}.pbs.perm", pc=range(1, 21), chr=CHR_VALUES),
        expand("output/01.qtltools_re/perm_fnf/pc{pc}/chr{chr}.fnf.perm", pc=range(1, 21), chr=CHR_VALUES),
        expand("output/01.qtltools_re/nominal_pbs_wasp/pc{pc}/chr{chr}.pbs.cis", pc=range(1, 21), chr=CHR_VALUES),
        expand("output/01.qtltools_re/nominal_fnf_wasp/pc{pc}/chr{chr}.fnf.cis", pc=range(1, 21), chr=CHR_VALUES),
        expand("output/01.qtltools_re/perm_pbs_wasp/pc{pc}/chr{chr}.pbs.perm", pc=range(1, 21), chr=CHR_VALUES),
        expand("output/01.qtltools_re/perm_fnf_wasp/pc{pc}/chr{chr}.fnf.perm", pc=range(1, 21), chr=CHR_VALUES),
        expand("output/01.qtltools_re/conditional_pbs/chr{chr}_pbs_condtional.txt",chr=CHR_VALUES),
        expand("output/01.qtltools_re/conditional_fnf/chr{chr}_fnf_condtional.txt",chr=CHR_VALUES)
        

rule create_vcf:
    output:
        vcf="output/geno/chQTL_samplesubset_101.vcf.gz"
    params:
        samples_to_exclude = config['vcf_sample_exluce_afterQC'],
        samtoolsVer = config['samtoolsVers'],
        raw_vcf = config['vcf']
    shell:
        """
        ml samtools/{params.samtoolsVer}
        mkdir -p output/geno
        bcftools view -s "^{params.samples_to_exclude}" {params.raw_vcf} -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}

        """

rule samplename_mapping:
    input:
        VCF=rules.create_vcf.output.vcf,
        PBS='/work/users/s/e/seyoun/CQTL_sQTL/output/CTL_samples.txt',
        FNF='/work/users/s/e/seyoun/CQTL_sQTL/output/FNF_samples.txt',
        OA='/work/users/s/e/seyoun/CQTL_sQTL/output/OA_samples.txt'
    output:
        pbs_matched="output/geno/pbs_geno/pbs_matched_afterQC.txt",
        fnf_matched="output/geno/fnf_geno/fnf_matched_afterQC.txt",
        oa_matched="output/geno/oa_geno/oa_matched_afterQC.txt"
    params:
        samtoolsVer=config['samtoolsVers'],
        VCF_OA=config['vcf_oa']

    log:
        out="output/logs/processVCF_sqtl.out",
        err="output/logs/processVCF_sqtl.err"

    shell:
        """
        ml samtools/{params.samtoolsVer}
        mkdir -p output/geno/pbs_geno
        mkdir -p output/geno/fnf_geno
        mkdir -p output/geno/oa_geno
       # Run PBS

        sh /work/users/s/e/seyoun/CQTL_sQTL/scripts/namechange.sh {input.VCF} {input.PBS}  {output.pbs_matched} 1> {log.out} 2> {log.err}

        # Run FNF

        sh /work/users/s/e/seyoun/CQTL_sQTL/scripts/namechange.sh {input.VCF} {input.FNF}  {output.fnf_matched}
        
         
        # Run OA

        sh /work/users/s/e/seyoun/CQTL_sQTL/scripts/namechange.sh {params.VCF_OA} {input.OA}  {output.oa_matched}
        
        """

rule re_processVCF:
    input:
        VCF=rules.create_vcf.output.vcf,
        VCFOA_in=config['vcf_oa'],
        vcf_pbs=rules.samplename_mapping.output.pbs_matched,
        vcf_fnf=rules.samplename_mapping.output.fnf_matched,
        vcf_oa=rules.samplename_mapping.output.oa_matched
    output:
        pbs_vcf="output/geno/pbs_geno/pbs_rename_101_wCHR.vcf.gz",
        fnf_vcf="output/geno/fnf_geno/fnf_rename_101_wCHR.vcf.gz",
        oa_vcf="output/geno/oa_geno/oa_rename_101_wCHR.vcf.gz",
        pbs_temp="output/geno/pbs_geno/pbs_rename_101.vcf.gz",
        fnf_temp="output/geno/fnf_geno/fnf_rename_101.vcf.gz",
        oa_temp="output/geno/oa_geno/oa_rename_101.vcf.gz",
        pbs_fi_vcf="output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz",
        fnf_fi_vcf='output/geno/fnf_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz',
        oa_fi_vcf='output/geno/oa_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz'

    params:
        sample_exclue=config['vcf_sample_exluce_afterQC'],
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

        sh /work/users/s/e/seyoun/CQTL_sQTL/scripts/vcfpreprocess.sh {input.VCF} {input.vcf_pbs} {output.pbs_temp} {output.pbs_vcf} {params.pbs_dir} 1> {log.pbsout} 2> {log.pbserr}
        
        sh /work/users/s/e/seyoun/CQTL_sQTL/scripts/vcfpreprocess.sh {input.VCF} {input.vcf_fnf} {output.fnf_temp} {output.fnf_vcf} {params.fnf_dir} 1> {log.fnfout} 2> {log.fnferr}

        sh /work/users/s/e/seyoun/CQTL_sQTL/scripts/vcfpreprocess.sh {input.VCF} {input.vcf_oa} {output.oa_temp} {output.oa_vcf} {params.oa_dir} 1> {log.oaout} 2> {log.oaerr}
        
        """

rule junclist:
    output:
        ctl_fnf = 'output/junc/ctl_fnf_junczipList.txt',
        ctl_fnf_wasp = 'output/junc_wasp/ctl_fnf_junczipList.txt'
    log:
        err= 'output/logs/junczip_text.err',
        err_wasp= 'output/logs/junczip_wasp_text.err'
    shell:
        """
        find /work/users/s/e/seyoun/CQTL_sQTL/output/junc -type f \( -name "*FNF*" -o -name "*CTL*" \) -name "*.junc*" -print > {output.ctl_fnf} 2> {log.err}
        find /work/users/s/e/seyoun/CQTL_sQTL/output/junc_wasp -type f \( -name "*FNF*" -o -name "*CTL*" \) -name "*.junc*" -print > {output.ctl_fnf_wasp} 2> {log.err_wasp}

        """

rule cluster_prepare:
    input:
        junclist=rules.junclist.output.ctl_fnf,
        junclist_wasp=rules.junclist.output.ctl_fnf_wasp
    output:
        cluster='output/gtex_cluster/ctl_fnf_perind.counts.filtered.gz',
        cluster_wasp='output/gtex_cluster_wasp/ctl_fnf_perind.counts.filtered.gz',
        bed= "output/gtex_cluster/ctl_fnf.leafcutter.bed.gz",
        bed_wasp="output/gtex_cluster_wasp/ctl_fnf.leafcutter.bed.gz",
        pca="output/gtex_cluster/ctl_fnf.leafcutter.PCs.txt",
        pca_wasp="output/gtex_cluster_wasp/ctl_fnf.leafcutter.PCs.txt"
        #cluster_dir ='output/gtex_cluster',
        #cluster_dir_wasp='output/gtex_cluster_wasp'
    log:
        err_genes="output/logs/gene.gtf.err",
        err_exon="output/logs/exon.gtf.err",
        err_clu="output/logs/cluster_gtx.err",
        err_clu_wasp="output/logs/cluster_gtx_wasp.err"
    params:
        pVers=config['pythonVers'],
        rVers=config['rVers'],
        gtexDir= config['gtex_dir'],
        gtf= config['gtf'],
        gtf_dir= config['genome_dir'],
        leafcutterDir=config['leafcutter']
    shell:
        """
        ml python/{params.pVers}
        ml r/{params.rVers}
        #mkdir -p output/gtex_cluster
        #mkdir -p output/gtex_cluster_wasp

        #python {params.gtexDir}/collapse_annotation.py {params.gtf} {params.gtf_dir}/gencode.v45.annotation.genes.gtf > {log.err_genes} 2>&1
        python {params.gtexDir}/make_exonlist.py > {log.err_exon} 2>&1
        python {params.gtexDir}/cluster_prepare_fastqtl.py {input.junclist} \
                {params.gtf_dir}/gencode.v45.annotation.genes.exons.txt.gz \
                {params.gtf_dir}/gencode.v45.annotation.genes.gtf \
                ctl_fnf \
                /work/users/s/e/seyoun/CQTL_sQTL/output/junc/sample_participant_lookup.txt \
                --min_clu_reads 100 \
                --max_intron_len 500000 \
                --num_pcs 20 \
                --leafcutter_dir {params.leafcutterDir} \
                -o output/gtex_cluster > {log.err_clu} 2>&1
        echo {output.cluster}       
        python {params.gtexDir}/cluster_prepare_fastqtl.py {input.junclist_wasp} \
                {params.gtf_dir}/gencode.v45.annotation.genes.exons.txt.gz \
                {params.gtf_dir}/gencode.v45.annotation.genes.gtf \
                ctl_fnf \
                /work/users/s/e/seyoun/CQTL_sQTL/output/junc_wasp/sample_participant_lookup.txt \
                --min_clu_reads 100 \
                --max_intron_len 500000 \
                --num_pcs 20 \
                --leafcutter_dir {params.leafcutterDir} \
                -o output/gtex_cluster_wasp > {log.err_clu_wasp} 2>&1
        echo {output.cluster_wasp}

        """


#rule leafcutter_qqnorm:
#    input:
#        count= rules.cluster_prepare.output.cluster,
#        count_wasp= rules.cluster_prepare.output.cluster_wasp
#    output:
#        qqnorm= '/output/gtex_cluster/ctl_fnf_perind.counts.filtered.gz_prepare.sh',
#        qqnorm_wasp= 'output/gtex_cluster_wasp/ctl_fnf_perind.counts.filtered.gz_prepare.sh'
#    log:
#        err1="output/logs/qqnorm_leafcutter.err",
#        err2="output/logs/qqnorm_wasp_leafcutter.err"
#    params:
#        pVers=config['pythonVers'],
#        rVers=config['rVers'],
#        leafcutterDir=config['leafcutter']
#    shell:
#        """
#        ml python/{params.pVers}
#        ml r/{params.rVers}
#        
#        python3 {params.leafcutterDir}scripts/prepare_phenotype_table.py {input.count} -p20 > {log.err1} 2>&1
#        #sh output/gtex_cluster/ctl_fnf_perind.counts.filtered.gz_prepare.sh
#        echo {output.qqnorm}
#        python3 {params.leafcutterDir}scripts/prepare_phenotype_table.py {input.count_wasp} -p20 > {log.err2} 2>&1
#        #sh output/gtex_cluster_wasp/ctl_fnf_perind.counts.filtered.gz_prepare.sh
#        #echo {output.qqnorm_wasp}
#
#        """


rule prep_qtltools:
    input:
        pca=rules.cluster_prepare.output.pca,
        pca_wasp=rules.cluster_prepare.output.pca_wasp
    output:
        cov= "output/gtex_cluster/qtltools_prep/covariates_PC1",
        cov_wasp="output/gtex_cluster_wasp/qtltools_prep/covariates_PC1"
    log:
        err1="output/logs/qtlprep_PC.err",
        err2="output/logs/qtlprep_wasp_PC.err"
    
    params:
        rVers=config['rVers'] 
    shell:
        """
        ml r/{params.rVers}
        Rscript scripts/sQTL_rscripts/bedfile_seperatebyCHR.R /work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster/ output/gtex_cluster/qtltools_prep {input.pca} > {log.err1} 2>&1
        sh output/gtex_cluster/qtltools_prep/tabix.sh

        Rscript scripts/sQTL_rscripts/bedfile_seperatebyCHR.R /work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster_wasp/ output/gtex_cluster_wasp/qtltools_prep {input.pca_wasp} wasp > {log.err2} 2>&1
        sh output/gtex_cluster_wasp/qtltools_prep/tabix.sh
        """


rule run_qtltools:
    input:
        vcf_pbs=rules.re_processVCF.output.pbs_fi_vcf,
        vcf_fnf=rules.re_processVCF.output.fnf_fi_vcf,
        cov="output/gtex_cluster/qtltools_prep/covariates_PC{pc}",
        bed="output/gtex_cluster/qtltools_prep/ctrvsfnf_qqnorm_chr{chr}.bed.gz"
    output:
        nominal_pbs="output/01.qtltools_re/nominal_pbs/pc{pc}/chr{chr}.pbs.cis",
        nominal_fnf="output/01.qtltools_re/nominal_fnf/pc{pc}/chr{chr}.fnf.cis",
        perm_pbs="output/01.qtltools_re/perm_pbs/pc{pc}/chr{chr}.pbs.perm",
        perm_fnf="output/01.qtltools_re/perm_fnf/pc{pc}/chr{chr}.fnf.perm"
    log:
        err1="output/logs/qtl_nominal_pbs_pc{pc}_chr{chr}.err",
        err2="output/logs/qtl_nominal_fnf_pc{pc}_chr{chr}.err",
        err3="output/logs/qtl_perm_pbs_pc{pc}_chr{chr}.err",
        err4="output/logs/qtl_perm_fnf_pc{pc}_chr{chr}.err"
    params:
        qtltoolsVer=config['qtltools']

    shell:
        """
        module load qtltools/{params.qtltoolsVer}

        mkdir -p output/01.qtltools_re/nominal_pbs/pc{wildcards.pc}
        mkdir -p output/01.qtltools_re/nominal_fnf/pc{wildcards.pc}
        mkdir -p output/01.qtltools_re/perm_pbs/pc{wildcards.pc}
        mkdir -p output/01.qtltools_re/perm_fnf/pc{wildcards.pc}

        QTLtools cis --vcf {input.vcf_pbs} --bed {input.bed} --cov {input.cov} --window 100000 --out {output.nominal_pbs} --nominal 1.0 --std-err > {log.err1} 2>&1
        QTLtools cis --vcf {input.vcf_fnf} --bed {input.bed} --cov {input.cov} --window 100000 --out {output.nominal_fnf} --nominal 1.0 --std-err > {log.err2} 2>&1
        QTLtools cis --vcf {input.vcf_pbs} --bed {input.bed} --cov {input.cov} --window 100000 --out {output.perm_pbs} --permute 1000 > {log.err3} 2>&1
        QTLtools cis --vcf {input.vcf_fnf} --bed {input.bed} --cov {input.cov} --window 100000 --out {output.perm_fnf} --permute 1000 > {log.err4} 2>&1
        """

rule run_qtltools_wasp:
    input:
        vcf_pbs=rules.re_processVCF.output.pbs_fi_vcf,
        vcf_fnf=rules.re_processVCF.output.fnf_fi_vcf,
        cov="output/gtex_cluster_wasp/qtltools_prep/covariates_PC{pc}",
        bed="output/gtex_cluster_wasp/qtltools_prep/ctrvsfnf_qqnorm_chr{chr}.bed.gz"
    output:
        nominal_pbs="output/01.qtltools_re/nominal_pbs_wasp/pc{pc}/chr{chr}.pbs.cis",
        nominal_fnf="output/01.qtltools_re/nominal_fnf_wasp/pc{pc}/chr{chr}.fnf.cis",
        perm_pbs="output/01.qtltools_re/perm_pbs_wasp/pc{pc}/chr{chr}.pbs.perm",
        perm_fnf="output/01.qtltools_re/perm_fnf_wasp/pc{pc}/chr{chr}.fnf.perm"
    log:
        err1="output/logs/qtl_nominal_pbs_wasp_pc{pc}_chr{chr}.err",
        err2="output/logs/qtl_nominal_fnf_wasp_pc{pc}_chr{chr}.err",
        err3="output/logs/qtl_perm_pbs_wasp_pc{pc}_chr{chr}.err",
        err4="output/logs/qtl_perm_fnf_wasp_pc{pc}_chr{chr}.err"
    params:
        qtltoolsVer=config['qtltools']

    shell:
        """
        module load qtltools/{params.qtltoolsVer}

        mkdir -p output/01.qtltools_re/nominal_pbs_wasp/pc{wildcards.pc}
        mkdir -p output/01.qtltools_re/nominal_fnf_wasp/pc{wildcards.pc}
        mkdir -p output/01.qtltools_re/perm_pbs_wasp/pc{wildcards.pc}
        mkdir -p output/01.qtltools_re/perm_fnf_wasp/pc{wildcards.pc}

        QTLtools cis --vcf {input.vcf_pbs} --bed {input.bed} --cov {input.cov} --window 100000 --out {output.nominal_pbs} --nominal 1.0 --std-err > {log.err1} 2>&1
        QTLtools cis --vcf {input.vcf_fnf} --bed {input.bed} --cov {input.cov} --window 100000 --out {output.nominal_fnf} --nominal 1.0 --std-err > {log.err2} 2>&1
        QTLtools cis --vcf {input.vcf_pbs} --bed {input.bed} --cov {input.cov} --window 100000 --out {output.perm_pbs} --permute 1000 --std-err > {log.err3} 2>&1
        QTLtools cis --vcf {input.vcf_fnf} --bed {input.bed} --cov {input.cov} --window 100000 --out {output.perm_fnf} --permute 1000 --std-err > {log.err4} 2>&1
        """

rule runFDR_sig:
    input:
        perm_pbs="output/01.qtltools_re/perm_pbs/pc5_allchr.pbs.perm",
        perm_fnf="output/01.qtltools_re/perm_fnf/pc4_allchr.fnf.perm"
    output:
        pbs_sig="output/01.qtltools_re/01.significant/pbs_0.05_pc5.significant.txt",
        pbs_thres="output/01.qtltools_re/01.significant/pbs_0.05_pc5.thresholds.txt",
        fnf_sig="output/01.qtltools_re/01.significant/fnf_0.05_pc4.significant.txt",
        fnf_thres="output/01.qtltools_re/01.significant/fnf_0.05_pc4.thresholds.txt"
    log:
        err1="output/logs/runFDR_pbs_0.05.err",
        err2="output/logs/runFDR_fnf_0.05.err"
    params:
        qtltoolsVer=config['qtltools'],
        rVers=config['rVers'],
        calc_R="scripts/sQTL_rscripts/runFDR_cis.R"
    shell:
        """
        module load qtltools/{params.qtltoolsVer}
        ml r/{params.rVers}

        mkdir -p output/01.qtltools_re/01.significant
        #gzip -c {input.perm_pbs}
        #gzip -c {input.perm_fnf}

        Rscript {params.calc_R} {input.perm_pbs} 0.05 output/01.qtltools_re/01.significant/pbs_0.05_pc5 > {log.err1} 2>&1
        Rscript {params.calc_R} {input.perm_fnf} 0.05 output/01.qtltools_re/01.significant/fnf_0.05_pc4 > {log.err2} 2>&1

        """

rule runConditional:
    input:
        thres_pbs=rules.runFDR_sig.output.pbs_thres,
        thres_fnf=rules.runFDR_sig.output.fnf_thres,
        vcf_pbs=rules.re_processVCF.output.pbs_fi_vcf,
        vcf_fnf=rules.re_processVCF.output.fnf_fi_vcf,
        cov_pbs="output/gtex_cluster/qtltools_prep/covariates_PC5",
        cov_fnf="output/gtex_cluster/qtltools_prep/covariates_PC4",
        bed="output/gtex_cluster/qtltools_prep/ctrvsfnf_qqnorm_chr{chr}.bed.gz"
    output:
        cond_pbs="output/01.qtltools_re/conditional_pbs/chr{chr}_pbs_condtional.txt",
        cond_fnf="output/01.qtltools_re/conditional_fnf/chr{chr}_fnf_condtional.txt"
    log:
        err1="output/logs/conditional_chr{chr}_pbs.err",
        err2="output/logs/conditional_chr{chr}_fnf.err"
    params:
        qtltoolsVer=config['qtltools']
    shell:
        """
        module load qtltools/{params.qtltoolsVer}

        mkdir -p output/01.qtltools_re/conditional_pbs
        mkdir -p output/01.qtltools_re/conditional_fnf;

 
        QTLtools cis --vcf {input.vcf_pbs} --bed {input.bed}  --cov {input.cov_pbs} --mapping {input.thres_pbs} --window 100000 --out {output.cond_pbs} > {log.err1} 2>&1
        QTLtools cis --vcf {input.vcf_fnf} --bed {input.bed}  --cov {input.cov_fnf} --mapping {input.thres_fnf} --window 100000 --out {output.cond_fnf} > {log.err2} 2>&1
       """ 
