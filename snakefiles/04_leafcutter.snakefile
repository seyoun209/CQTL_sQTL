#!/usr/bin/env python3
# -*- coding: utf-8 -*-

rule all:
    input:
        'output/clu_fnf/PBSvsFNF.Rdata',
        'output/clu_oa/PBSvsOA.Rdata',
        'output/clu_fnf/ctlvsfnf_perind.counts.gz.PCs',
        'output/clu_oa/ctlvsoa_perind.counts.gz.PCs'

rule junclist:
    output:
        ctl_fnf = 'output/junc/ctl_fnf_juncList.txt',
        ctl_oa = 'output/junc/ctl_oa_juncList.txt'
    log:
        err = 'output/logs/junc_text.err'
    shell:
        """
        find output/junc -type f \\( -name "*FNF*" -o -name "*CTL*" \\) -name "*.junc*" -print > {output.ctl_fnf} 2> {log.err}
        find output/junc -type f \\( -name "*OA*" -o -name "*CTL*" \\) -name "*.junc*" -print > {output.ctl_oa}
        """

rule cluster:
    input:
        set1 = rules.junclist.output.ctl_fnf,
        set2 = rules.junclist.output.ctl_oa
    output:
        counts1 = 'output/clu_fnf/ctlvsfnf_perind_numers.counts.gz',
        counts2 = 'output/clu_oa/ctlvsoa_perind_numers.counts.gz',
        ct1 = 'output/clu_fnf/ctlvsfnf_perind.counts.gz',
        ct2 = 'output/clu_oa/ctlvsoa_perind.counts.gz'
    params:
        pVers = config['pythonVers'],
        leafcutterVer = config['leafcutter']
    shell:
        """
        module load python/{params.pVers}
        mkdir -p output/clu_fnf output/clu_oa

        # FNF vs CTL
        python {params.leafcutterVer}/clustering/leafcutter_cluster.py -j {input.set1} -r output/clu_fnf
        gzip -c output/clu_fnf/leafcutter_perind_numers.counts > {output.counts1}
        gzip -c output/clu_fnf/leafcutter_perind.counts > {output.ct1}

        # OA vs CTL
        python {params.leafcutterVer}/clustering/leafcutter_cluster.py -j {input.set2} -r output/clu_oa
        gzip -c output/clu_oa/leafcutter_perind_numers.counts > {output.counts2}
        gzip -c output/clu_oa/leafcutter_perind.counts > {output.ct2}
        """

rule PCs:
    input:
        ct1 = rules.cluster.output.ct1,
        ct2 = rules.cluster.output.ct2
    output:
        pcs1 = 'output/clu_fnf/ctlvsfnf_perind.counts.gz.PCs',
        pcs2 = 'output/clu_oa/ctlvsoa_perind.counts.gz.PCs'
    params:
        rVers = config['rVers'],
        script = 'scripts/calculate_PCs.R'
    shell:
        """
        module load R/{params.rVers}
        Rscript {params.script} {input.ct1} {output.pcs1}
        Rscript {params.script} {input.ct2} {output.pcs2}
        """

rule leafcutter_ds:
    input:
        pcs1 = rules.PCs.output.pcs1,
        pcs2 = rules.PCs.output.pcs2
    output:
        rdata1 = 'output/clu_fnf/PBSvsFNF.Rdata',
        rdata2 = 'output/clu_oa/PBSvsOA.Rdata'
    params:
        rVers = config['rVers'],
        script = 'scripts/leafcutter_ds.R'
    shell:
        """
        module load R/{params.rVers}
        Rscript {params.script} {input.pcs1} {output.rdata1}
        Rscript {params.script} {input.pcs2} {output.rdata2}
        """
