#!/usr/bin/env python3

rule categorize_samples:
    output:
        ctl="output/CTL_samples.txt",
        fnf="output/FNF_samples.txt",
        oa="output/OA_samples.txt"
    shell:
        """
        sh /work/users/s/e/seyoun/CQTL_sQTL/scripts/categorize_samples.sh
        """

rule create_vcf:
    output:
        vcf="output/geno/chQTL_samplesubset_103.vcf.gz"
    params:
        samples_to_exclude = config['vcf_sample_exclude'],
        samtoolsVer = config['samtoolsVers'],
        raw_vcf = config['vcf']
    shell:
        """
        ml samtools/{params.samtoolsVer}
        mkdir -p output/geno
        bcftools view -s "^{params.samples_to_exclude}" {params.raw_vcf} -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}

        """

rule processVCF:
    input:
        VCF=rules.create_vcf.output.vcf,
        PBS=rules.categorize_samples.output.ctl,
        FNF=rules.categorize_samples.output.fnf,
        OA=rules.categorize_samples.output.oa
    output:
        pbs_matched="output/geno/pbs_geno/pbs_matched.txt",
        fnf_matched="output/geno/fnf_geno/fnf_matched.txt",
        oa_matched="output/geno/oa_geno/oa_matched.txt"
    params:
        samtoolsVer=config['samtoolsVers'],
        VCF_OA=config['vcf_oa']
    
    log:
        out="output/logs/processVCF.out",
        err="output/logs/processVCF.err"

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
        oa_temp="output/geno/oa_geno/oa.rename.vcf.gz"

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

        
        """


