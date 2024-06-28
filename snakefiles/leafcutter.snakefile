#!/usr/bin/env python3
# -*- coding: utf-8 -*-

rule all:
    input:
        'output/clu_fnf/PBSvsFNF.Rdata',
        'output/clu_oa/PBSvsOA.Rdata',
        'output/clu_fnf_wasp/PBSvsFNF_wasp.Rdata',
        'output/clu_oa_wasp/PBSvsOA_wasp.Rdata',
        [expand("output/{dir_nm}/ctlvsfnf_perind.counts.gz.PCs",dir_nm=['clu_fnf','clu_fnf_wasp'])]


rule junclist:
    output:
        ctl_fnf = 'output/junc/ctl_fnf_juncList.txt',
        ctl_oa = 'output/junc/ctl_oa_juncList.txt',
        ctl_fnf_wasp = 'output/junc_wasp/ctl_fnf_juncList.txt',
        ctl_oa_wasp = 'output/junc_wasp/ctl_oa_juncList.txt'
    log:
        err= 'output/logs/junc_text.err',
        err_wasp= 'output/logs/junc_wasp_text.err'
    shell:
        """
        find /work/users/s/e/seyoun/CQTL_sQTL/output/junc -type f \( -name "*FNF*" -o -name "*CTL*" \) -name "*.junc*" -print > {output.ctl_fnf} 2> {log.err}
        find /work/users/s/e/seyoun/CQTL_sQTL/output/junc -type f \( -name "*OA*" -o -name "*CTL*" \) -name "*.junc*" -print > {output.ctl_oa}
        find /work/users/s/e/seyoun/CQTL_sQTL/output/junc_wasp -type f \( -name "*FNF*" -o -name "*CTL*" \) -name "*.junc*" -print > {output.ctl_fnf_wasp} 2> {log.err_wasp}
        find /work/users/s/e/seyoun/CQTL_sQTL/output/junc_wasp -type f \( -name "*OA*" -o -name "*CTL*" \) -name "*.junc*" -print > {output.ctl_oa_wasp}
        """

rule cluster:
    input:
        set1=rules.junclist.output.ctl_fnf,
        set2=rules.junclist.output.ctl_oa,
        set3=rules.junclist.output.ctl_fnf_wasp,
        set4=rules.junclist.output.ctl_oa_wasp
    output:
        counts1='output/clu_fnf/ctlvsfnf_perind_numers.counts.gz',
        counts2='output/clu_oa/ctlvsoa_perind_numers.counts.gz',
        counts3= 'output/clu_fnf_wasp/ctlvsfnf_perind_numers.counts.gz',
        counts4= 'output/clu_oa_wasp/ctlvsoa_perind_numers.counts.gz',
        ct1='output/clu_fnf/ctlvsfnf_perind.counts.gz',
        ct2='output/clu_oa/ctlvsoa_perind.counts.gz',
        ct3='output/clu_fnf_wasp/ctlvsfnf_perind.counts.gz',
        ct4='output/clu_oa_wasp/ctlvsoa_perind.counts.gz'

    params:
        pVers=config['pythonVers'],
        rVers=config['rVers'],
        leafcutterVer=config['leafcutter']
    log:
        err1='output/logs/leafcutter_ctl_fnf.out',
        err2='output/logs/leafcutter_ctl_oa.out',
        err3='output/logs/leafcutter_wasp_ctl_fnf.out',
        err4='output/logs/leafcutter_wasp_ctl_oa.out'
    shell:
        """
        ml python/{params.pVers}
        ml r/{params.rVers}
        mkdir -p output/clu_fnf
        mkdir -p output/clu_oa
        mkdir -p output/clu_fnf_wasp
        mkdir -p output/clu_oa_wasp
        mkdir -p output/clu_fnf/perind
        mkdir -p output/clu_oa/perind
        mkdir -p output/clu_fnf_wasp/perind
        mkdir -p output/clu_oa_wasp/perind

        python {params.leafcutterVer}clustering/leafcutter_cluster_regtools.py -j {input.set1} -m 100 -o output/clu_fnf/ctlvsfnf -l 500000 > {log.err1} 2>&1
        mv *sorted.gz* output/clu_fnf/perind 
        python {params.leafcutterVer}clustering/leafcutter_cluster_regtools.py -j {input.set2} -m 100 -o output/clu_oa/ctlvsoa -l 500000 > {log.err2} 2>&1
        mv *.gz* output/clu_oa/perind
        python {params.leafcutterVer}clustering/leafcutter_cluster_regtools.py -j {input.set3} -m 100 -o output/clu_fnf_wasp/ctlvsfnf -l 500000  > {log.err3} 2>&1
        mv *.gz* output/clu_fnf_wasp/perind
        python {params.leafcutterVer}clustering/leafcutter_cluster_regtools.py -j {input.set4} -m 100 -o output/clu_oa_wasp/ctlvsoa -l 500000  > {log.err4} 2>&1
        mv *.gz* output/clu_oa_wasp/perind

        """

rule diff_leafcutter:
    input:
        ct1 = rules.cluster.output.counts1,
        ct2 = rules.cluster.output.counts2,
        ct3 = rules.cluster.output.counts3,
        ct4 = rules.cluster.output.counts4
    output:
        ds1='output/clu_fnf/ctlvsfnf_ds_cluster_significance.txt',
        ds2='output/clu_oa/ctlvsoa_ds_cluster_significance.txt',
        ds3='output/clu_fnf_wasp/ctlvsfnf_wasp_ds_cluster_significance.txt',
        ds4='output/clu_oa_wasp/ctlvsoa_was_ds_cluster_significance.txt',
        es1='output/clu_fnf/ctlvsfnf_ds_effect_sizes.txt',
        es2='output/clu_oa/ctlvsoa_ds_effect_sizes.txt',
        es3='output/clu_fnf_wasp/ctlvsfnf_wasp_ds_effect_sizes.txt',
        es4='output/clu_oa_wasp/ctlvsoa_was_ds_effect_sizes.txt'
    log:
        err_txt1='output/logs/leafcutter_group_fnf.err',
        err_txt2='output/logs/leafcutter_group_oa.err',
        err_txt3='output/logs/leafcutter_group_fnf_wasP.err',
        err_txt4='output/logs/leafcutter_group_oa_wasp.err',
        err_ds1='output/logs/leafcutter_ds_ctrvsfnf.err',
        err_ds2='output/logs/leafcutter_ds_ctrvsoa.err',
        err_ds3='output/logs/leafcutter_ds_wasp_ctrvsfnf.err',
        err_ds4='output/logs/leafcutter_ds_wasp_ctrvsoa.err'


    params:
        pVers=config['pythonVers'],
        rVers=config['rVers'],
        leafcutterVer=config['leafcutter']
    threads: 6
    shell:
        """
        ml r/{params.rVers}
        ml python/{params.pVers}
        Rscript scripts/process_samples.R ./aligned_samplesheet.txt "./donor_samples.txt" "./rna_extraction.txt" "CTL FNF" > {log.err_txt1} 2>&1
        Rscript scripts/process_samples.R ./aligned_samplesheet.txt "./donor_samples.txt" "./rna_extraction.txt" "CTL OA"  > {log.err_txt2} 2>&1
       Rscript scripts/process_samples.R ./aligned_samplesheet.txt "./donor_samples.txt" "./rna_extraction.txt" "CTL FNF" wasp > {log.err_txt3} 2>&1
        Rscript scripts/process_samples.R ./aligned_samplesheet.txt "./donor_samples.txt" "./rna_extraction.txt" "CTL OA" wasp > {log.err_txt4} 2>&1
 


        Rscript {params.leafcutterVer}scripts/leafcutter_ds.R --num_threads 6 {input.ct1} output/CTLvsFNF_group.txt -i 10 -g 10 -e {params.leafcutterVer}gencode_hg38_exon_v2.txt.gz -o output/clu_fnf/ctlvsfnf_ds > {log.err_ds1} 2>&1
        Rscript {params.leafcutterVer}scripts/leafcutter_ds.R --num_threads 6 {input.ct2} output/CTLvsOA_group.txt -e {params.leafcutterVer}gencode_hg38_exon_v2.txt.gz -o output/clu_oa/ctlvsoa_ds > {log.err_ds2} 2>&1
        Rscript {params.leafcutterVer}scripts/leafcutter_ds.R --num_threads 6 {input.ct3} output/wasp_CTLvsFNF_group.txt -i 10 -g 10 -e {params.leafcutterVer}gencode_hg38_exon_v2.txt.gz -o output/clu_fnf_wasp/ctlvsfnf_wasp_ds > {log.err_ds3} 2>&1
        Rscript {params.leafcutterVer}scripts/leafcutter_ds.R --num_threads 6 {input.ct4} output/wasp_CTLvsOA_group.txt -e {params.leafcutterVer}gencode_hg38_exon_v2.txt.gz -o output/clu_oa_wasp/ctlvsoa_was_ds > {log.err_ds4} 2>&1


        """

rule leafviz:
    input:
        ds1=rules.diff_leafcutter.output.ds1,
        es1=rules.diff_leafcutter.output.es1,
        ds2=rules.diff_leafcutter.output.ds2,
        es2=rules.diff_leafcutter.output.es2,
        ds3=rules.diff_leafcutter.output.ds3,
        es3=rules.diff_leafcutter.output.es3,
        ds4=rules.diff_leafcutter.output.ds4,
        es4=rules.diff_leafcutter.output.es4,
        ct1 = rules.cluster.output.counts1,
        ct2 = rules.cluster.output.counts2,
        ct3 = rules.cluster.output.counts3,
        ct4 = rules.cluster.output.counts4

    output:
        rdata1 = 'output/clu_fnf/PBSvsFNF.Rdata',
        rdata2 = 'output/clu_oa/PBSvsOA.Rdata',
        rdata3 = 'output/clu_fnf_wasp/PBSvsFNF_wasp.Rdata',
        rdata4 = 'output/clu_oa_wasp/PBSvsOA_wasp.Rdata'
    log:
        err1='output/logs/leafviz_ctr_fnf.err',
        err2='output/logs/leafviz_ctr_oa.err',
        err3='output/logs/leafviz_ctr_fnf_wasp.err',
        err4='output/logs/leafviz_ctr_oa_wasp.err'
    
    params:
        pVers=config['pythonVers'],
        rVers=config['rVers'],
        leafcutterVer=config['leafcutter']
    shell:
        """
        ml r/{params.rVers}
        ml python/{params.pVers}

        Rscript {params.leafcutterVer}leafviz/prepare_results.R \
                -m output/CTLvsFNF_group.txt \
                --code leafcutter {input.ct1} {input.ds1} {input.es1} \
                {params.leafcutterVer}gencode_hg38 \
                -o {output.rdata1} > {log.err1} 2>&1
        Rscript {params.leafcutterVer}leafviz/prepare_results.R \
                -m output/CTLvsOA_group.txt \
                --code leafcutter {input.ct2} {input.ds2} {input.es2}\
                {params.leafcutterVer}gencode_hg38 \
                -o {output.rdata2} > {log.err2} 2>&1

        Rscript {params.leafcutterVer}leafviz/prepare_results.R -m output/wasp_CTLvsFNF_group.txt --code leafcutter {input.ct3} {input.ds3} {input.es3} {params.leafcutterVer}gencode_hg38 -o {output.rdata3} > {log.err3} 2>&1
        Rscript {params.leafcutterVer}leafviz/prepare_results.R -m output/wasp_CTLvsOA_group.txt --code leafcutter {input.ct4} {input.ds4} {input.es4} {params.leafcutterVer}gencode_hg38 -o {output.rdata4} > {log.err4} 2>&1

        """

rule pca:
    input:
        perind = rules.cluster.output.ct1,
        perind_wasp = rules.cluster.output.ct3
    output:
        qqnorm = 'output/clu_fnf/ctlvsfnf_perind.counts.gz.PCs',
        pca_wasp = 'output/clu_fnf_wasp/ctlvsfnf_perind.counts.gz.PCs'
    log:
        err1='output/logs/pc_calc_fnf.err',
        err2='output/logs/pc_calc_fnf_wasp.err'
    params:
        pVers=config['pythonVers'],
        rVers=config['rVers'],
        leafcutterVer=config['leafcutter']
    shell:
        """
        ml r/{params.rVers}
        ml python/{params.pVers}

        python {params.leafcutterVer}scripts/prepare_phenotype_table.py {input.perind} -p 20 > {log.err1} 2>&1
        python {params.leafcutterVer}scripts/prepare_phenotype_table.py {input.perind_wasp} -p 20 > {log.err2} 2>&1

        """

