#!/usr/bin/env python3

CHR_VALUES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10","11", "12", "13", "14", "15", "16", "17", "18", "19","20", "21", "22"]

rule all:
    input:
        "output/geno/oa_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz",
        "output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz",
        "output/junc/sample_participant_lookup_pbs_oa.txt",
        cluster_prepare=expand("output/gtex_cluster_oa/pbs_oa{ext}",ext=['_perind.counts.filtered.gz','.leafcutter.bed.gz','.leafcutter.PCs.txt','_perind.counts.filtered.gz_prepare.sh']),
        covariates=expand("output/gtex_cluster_oa/qtltools_prep/covariates_PC{pc}", pc=range(1,21)),
        bed_files=expand("output/gtex_cluster_oa/qtltools_prep/pbs_oa_qqnorm_chr{chr}.bed.gz", chr=CHR_VALUES),
        bed_indexes=expand("output/gtex_cluster_oa/qtltools_prep/pbs_oa_qqnorm_chr{chr}.bed.gz.tbi", chr=CHR_VALUES),
        # QTL outputs - split into 4 separate expand calls for correct file extensions
        nominal_pbs=expand("output/03.qtltools_re_oa/nominal_pbs/pc{pc}/chr{chr}.pbs.cis", 
                          pc=range(1,21), chr=CHR_VALUES),
        nominal_oa=expand("output/03.qtltools_re_oa/nominal_oa/pc{pc}/chr{chr}.oa.cis", 
                         pc=range(1,21), chr=CHR_VALUES),
        perm_pbs=expand("output/03.qtltools_re_oa/perm_pbs/pc{pc}/chr{chr}.pbs.perm", 
                       pc=range(1,21), chr=CHR_VALUES),
        perm_oa=expand("output/03.qtltools_re_oa/perm_oa/pc{pc}/chr{chr}.oa.perm", 
                      pc=range(1,21), chr=CHR_VALUES),
        
        # Combined results - also separate expand calls for correct extensions
        nominal_pbs_combined=expand("output/03.qtltools_re_oa/nominal_pbs/pc{pc}_allchr.pbs.cis", 
                                  pc=range(1,21)),
        nominal_oa_combined=expand("output/03.qtltools_re_oa/nominal_oa/pc{pc}_allchr.oa.cis", 
                                 pc=range(1,21)),
        perm_pbs_combined=expand("output/03.qtltools_re_oa/perm_pbs/pc{pc}_allchr.pbs.perm", 
                               pc=range(1,21)),
        perm_oa_combined=expand("output/03.qtltools_re_oa/perm_oa/pc{pc}_allchr.oa.perm", 
                              pc=range(1,21)),
        # FDR results for PC6
        fdr_results=[
            "output/03.qtltools_re_oa/01.significant/pbs_0.05_pc6.significant.txt",
            "output/03.qtltools_re_oa/01.significant/pbs_0.05_pc6.thresholds.txt",
            "output/03.qtltools_re_oa/01.significant/oa_0.05_pc6.significant.txt",
            "output/03.qtltools_re_oa/01.significant/oa_0.05_pc6.thresholds.txt"
        ]

# Create PBS vs OA junclist
rule pbs_oa_junclist:
    output:
        pbs_oa = 'output/junc/pbs_oa_junczipList.txt'
    log:
        err= 'output/logs/pbs_oa_junczip_text.err'
    shell:
        """
        find /work/users/s/e/seyoun/CQTL_sQTL/output/junc -type f \( -name "*CTL*" -o -name "*OA*" \) -name "*.junc*" -print > {output.pbs_oa} 2> {log.err}
        """

# Create sample-participant lookup file for PBS and OA samples
rule create_pbs_oa_lookup:
    output:
        lookup = 'output/junc/sample_participant_lookup_pbs_oa.txt'
    shell:
        """
        # Find all PBS and OA sample files and extract the sample IDs
        for file in $(find /work/users/s/e/seyoun/CQTL_sQTL/output/junc -type f \( -name "*CTL*" -o -name "*OA*" \) -name "*.junc*"); do
            # Extract just the filename without path and extension
            SAMPLE=$(basename "$file" .junc.gz)
            echo -e "$SAMPLE\t$SAMPLE" >> {output.lookup}
        done
        
        # Remove duplicates while preserving order
        sort -u {output.lookup} -o {output.lookup}
        
        # Make sure header is at the top
        sed -i '1i Sample_ID\tParticipant_ID' {output.lookup}
        """

# Prepare PBS vs OA clusters
rule pbs_oa_cluster_prepare:
    input:
        junclist=rules.pbs_oa_junclist.output.pbs_oa,
        lookup=rules.create_pbs_oa_lookup.output.lookup
    output:
        cluster='output/gtex_cluster_oa/pbs_oa_perind.counts.filtered.gz',
        bed="output/gtex_cluster_oa/pbs_oa.leafcutter.bed.gz",
        pca="output/gtex_cluster_oa/pbs_oa.leafcutter.PCs.txt",
        qqnorm_script="output/gtex_cluster_oa/pbs_oa_perind.counts.filtered.gz_prepare.sh"
    log:
        err_clu="output/logs/cluster_oa_gtx.err"
    params:
        pVers=config['pythonVers'],
        rVers=config['rVers'],
        gtexDir= config['gtex_dir'],
        gtf_dir= config['genome_dir'],
        leafcutterDir=config['leafcutter']
    shell:
        """
        ml python/{params.pVers}
        ml r/{params.rVers}
        mkdir -p output/gtex_cluster_oa

        # Run the cluster_prepare_fastqtl.py script
        python {params.gtexDir}/cluster_prepare_fastqtl.py {input.junclist} \
                {params.gtf_dir}/gencode.v45.annotation.genes.exons.txt.gz \
                {params.gtf_dir}/gencode.v45.annotation.genes.gtf \
                pbs_oa \
                {input.lookup} \
                --min_clu_reads 100 \
                --max_intron_len 500000 \
                --num_pcs 20 \
                --leafcutter_dir {params.leafcutterDir} \
                -o output/gtex_cluster_oa > {log.err_clu} 2>&1

        # Generate the qqnorm script
        python {params.leafcutterDir}/scripts/prepare_phenotype_table.py \
                output/gtex_cluster_oa/pbs_oa_perind.counts.filtered.gz -p 20

        # Run the qqnorm script to create normalized files
        bash output/gtex_cluster_oa/pbs_oa_perind.counts.filtered.gz_prepare.sh

        # Create a backup of the qqnorm files
        mkdir -p output/gtex_cluster_oa/qqnorm_backup
        cp output/gtex_cluster_oa/pbs_oa_perind.counts.filtered.gz.qqnorm_chr* output/gtex_cluster_oa/qqnorm_backup/ || true

        echo {output.cluster}
        """

# Prepare PBS vs OA QTLtools inputs
rule prep_pbs_oa_qtltools:
    input:
        pca=rules.pbs_oa_cluster_prepare.output.pca
    output:
        covs=expand("output/gtex_cluster_oa/qtltools_prep/covariates_PC{pc}", pc=range(1,21)),
        bed_files=expand("output/gtex_cluster_oa/qtltools_prep/pbs_oa_qqnorm_chr{chr}.bed.gz", chr=CHR_VALUES),
        bed_indexes=expand("output/gtex_cluster_oa/qtltools_prep/pbs_oa_qqnorm_chr{chr}.bed.gz.tbi", chr=CHR_VALUES),
        tabix_script="output/gtex_cluster_oa/qtltools_prep/tabix.sh"
    log:
        err="output/logs/qtlprep_pbs_oa.err"
    params:
        rVers=config['rVers']
    shell:
        """
        ml r/{params.rVers}
        mkdir -p output/gtex_cluster_oa/qtltools_prep

        Rscript scripts/revision/sQTL_OA/bedfile_seperatedbyCHR.R /work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster_oa/ /work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster_oa/qtltools_prep {input.pca} pbs_oa > {log.err} 2>&1
        
        sh {output.tabix_script}
        
        # Touch all output files to ensure timestamps are updated
        for f in {output.covs}; do
            touch -c $f
        done

        """

# Run QTLtools for PBS vs OA
rule run_qtltools_pbs_oa:
    input:
        vcf_pbs="output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz",
        vcf_oa="output/geno/oa_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz",
        cov="output/gtex_cluster_oa/qtltools_prep/covariates_PC{pc}",
        bed="output/gtex_cluster_oa/qtltools_prep/pbs_oa_qqnorm_chr{chr}.bed.gz",
        bed_index="output/gtex_cluster_oa/qtltools_prep/pbs_oa_qqnorm_chr{chr}.bed.gz.tbi",
        # Make sure all covariates are created before running QTLtools
        all_covs=rules.prep_pbs_oa_qtltools.output.covs
    output:
        nominal_pbs="output/03.qtltools_re_oa/nominal_pbs/pc{pc}/chr{chr}.pbs.cis",
        nominal_oa="output/03.qtltools_re_oa/nominal_oa/pc{pc}/chr{chr}.oa.cis",
        perm_pbs="output/03.qtltools_re_oa/perm_pbs/pc{pc}/chr{chr}.pbs.perm",
        perm_oa="output/03.qtltools_re_oa/perm_oa/pc{pc}/chr{chr}.oa.perm"
    log:
        err1="output/logs/qtl_nominal_oa_forPBS_pc{pc}_chr{chr}.err",
        err2="output/logs/qtl_nominal_oa_forOA_pc{pc}_chr{chr}.err",
        err3="output/logs/qtl_perm_oa_forPBS_pc{pc}_chr{chr}.err",
        err4="output/logs/qtl_perm_oa_forOA_pc{pc}_chr{chr}.err"
    params:
        qtltoolsVer=config['qtltools']
    shell:
        """
        module load qtltools/{params.qtltoolsVer}

        mkdir -p output/03.qtltools_re_oa/nominal_pbs/pc{wildcards.pc}
        mkdir -p output/03.qtltools_re_oa/nominal_oa/pc{wildcards.pc}
        mkdir -p output/03.qtltools_re_oa/perm_pbs/pc{wildcards.pc}
        mkdir -p output/03.qtltools_re_oa/perm_oa/pc{wildcards.pc}

        # Run QTLtools for PBS samples
        QTLtools cis --vcf {input.vcf_pbs} --bed {input.bed} --cov {input.cov} --window 100000 --out {output.nominal_pbs} --nominal 1.0 --std-err > {log.err1} 2>&1
        QTLtools cis --vcf {input.vcf_pbs} --bed {input.bed} --cov {input.cov} --window 100000 --out {output.perm_pbs} --permute 1000 > {log.err3} 2>&1
        
        # Run QTLtools for OA samples
        QTLtools cis --vcf {input.vcf_oa} --bed {input.bed} --cov {input.cov} --window 100000 --out {output.nominal_oa} --nominal 1.0 --std-err > {log.err2} 2>&1
        QTLtools cis --vcf {input.vcf_oa} --bed {input.bed} --cov {input.cov} --window 100000 --out {output.perm_oa} --permute 1000 > {log.err4} 2>&1
        """


# Rule to combine QTL results across chromosomes
rule combine_qtl_results:
    input:
        nominal_pbs = expand("output/03.qtltools_re_oa/nominal_pbs/pc{{pc}}/chr{chr}.pbs.cis", chr=CHR_VALUES),
        nominal_oa = expand("output/03.qtltools_re_oa/nominal_oa/pc{{pc}}/chr{chr}.oa.cis", chr=CHR_VALUES),
        perm_pbs = expand("output/03.qtltools_re_oa/perm_pbs/pc{{pc}}/chr{chr}.pbs.perm", chr=CHR_VALUES),
        perm_oa = expand("output/03.qtltools_re_oa/perm_oa/pc{{pc}}/chr{chr}.oa.perm", chr=CHR_VALUES)
    output:
        nominal_pbs_all = "output/03.qtltools_re_oa/nominal_pbs/pc{pc}_allchr.pbs.cis",
        nominal_oa_all = "output/03.qtltools_re_oa/nominal_oa/pc{pc}_allchr.oa.cis",
        perm_pbs_all = "output/03.qtltools_re_oa/perm_pbs/pc{pc}_allchr.pbs.perm",
        perm_oa_all = "output/03.qtltools_re_oa/perm_oa/pc{pc}_allchr.oa.perm"
    shell:
        """
        # Create directories if they don't exist
        mkdir -p output/03.qtltools_re_oa/nominal_pbs/
        mkdir -p output/03.qtltools_re_oa/nominal_oa/
        mkdir -p output/03.qtltools_re_oa/perm_pbs/
        mkdir -p output/03.qtltools_re_oa/perm_oa/
        
        # Combine PBS nominal results
        cat {input.nominal_pbs} > {output.nominal_pbs_all}
        
        # Combine OA nominal results
        cat {input.nominal_oa} > {output.nominal_oa_all}
        
        # Combine PBS permutation results
        cat {input.perm_pbs} > {output.perm_pbs_all}
        
        # Combine OA permutation results
        cat {input.perm_oa} > {output.perm_oa_all}
        """

rule runFDR_sig_pbs_oa:
    input:
        perm_pbs="output/03.qtltools_re_oa/perm_pbs/pc6_allchr.pbs.perm",
        perm_oa="output/03.qtltools_re_oa/perm_oa/pc6_allchr.oa.perm"
    output:
        pbs_sig="output/03.qtltools_re_oa/01.significant/pbs_0.05_pc6.significant.txt",
        pbs_thres="output/03.qtltools_re_oa/01.significant/pbs_0.05_pc6.thresholds.txt",
        oa_sig="output/03.qtltools_re_oa/01.significant/oa_0.05_pc6.significant.txt",
        oa_thres="output/03.qtltools_re_oa/01.significant/oa_0.05_pc6.thresholds.txt"
    log:
        err_pbs="output/logs/runFDR_pbs_pc6_0.05.err",
        err_oa="output/logs/runFDR_oa_pc6_0.05.err"
    params:
        qtltoolsVer=config['qtltools'],
        calc_R="scripts/sQTL_rscripts/runFDR_cis.R"
    shell:
        """
        module load qtltools/{params.qtltoolsVer}
        ml r/4.4.0

        mkdir -p output/03.qtltools_re_oa/01.significant

        # Run FDR for PBS
        Rscript {params.calc_R} {input.perm_pbs} 0.05 output/03.qtltools_re_oa/01.significant/pbs_0.05_pc6 > {log.err_pbs} 2>&1

        # Run FDR for OA
        Rscript {params.calc_R} {input.perm_oa} 0.05 output/03.qtltools_re_oa/01.significant/oa_0.05_pc6 > {log.err_oa} 2>&1
        """

## Run conditional analysis for PBS vs OA
#rule runConditional_pbs_oa:
#    input:
#        thres_pbs=rules.runFDR_sig_pbs_oa.output.pbs_thres,
#        thres_oa=rules.runFDR_sig_pbs_oa.output.oa_thres,
#        vcf_pbs="output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz",
#        vcf_oa="output/geno/oa_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz",
#        cov="output/gtex_cluster_oa/qtltools_prep/covariates_PC5",
#        bed="output/gtex_cluster_oa/qtltools_prep/pbs_oa_qqnorm_chr{chr}.bed.gz"
#    output:
#        cond_pbs="output/03.qtltools_re_oa/conditional_pbs_vs_oa/chr{chr}_pbs_condtional.txt",
#        cond_oa="output/03.qtltools_re_oa/conditional_oa/chr{chr}_oa_condtional.txt"
#    log:
#        err_pbs="output/logs/conditional_chr{chr}_pbs_vs_oa.err",
#        err_oa="output/logs/conditional_chr{chr}_oa.err"
#    params:
#        qtltoolsVer=config['qtltools']
#    shell:
#        """
#        module load qtltools/{params.qtltoolsVer}
#
#        mkdir -p output/03.qtltools_re_oa/conditional_pbs_vs_oa
#        mkdir -p output/03.qtltools_re_oa/conditional_oa
#
#        # Run conditional analysis for PBS
#        QTLtools cis --vcf {input.vcf_pbs} --bed {input.bed} --cov {input.cov} --mapping {input.thres_pbs} --window 100000 --out {output.cond_pbs} > {log.err_pbs} 2>&1
#
#        # Run conditional analysis for OA
#        QTLtools cis --vcf {input.vcf_oa} --bed {input.bed} --cov {input.cov} --mapping {input.thres_oa} --window 100000 --out {output.cond_oa} > {log.err_oa} 2>&1
#        """
