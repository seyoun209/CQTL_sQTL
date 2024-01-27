#!/usr/bin/env python3

# Rule to declare all targets
rule all:
    input:
        [expand("output/{condition}_samples.txt",condition=['CTL','FNF','OA'])],
        "output/geno/chQTL_samplesubset_103.vcf.gz",
        [expand("output/geno/{condition}_geno/{condition}_matched.txt",condition=['pbs','fnf','oa'])],
        [expand("output/geno/{condition}_geno/{condition}.rename_wCHR.vcf.gz",condition=['pbs','fnf','oa'])]


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
        samtoolsVer=config['samtoolsVers']
    
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

        sh /work/users/s/e/seyoun/CQTL_sQTL/scripts/namechange.sh {input.VCF} {input.PBS} > {output.pbs_matched} 1> {log.out} 2> {log.err}

        # Run FNF

        sh /work/users/s/e/seyoun/CQTL_sQTL/scripts/namechange.sh {input.VCF} {input.FNF} > {output.fnf_matched}

        # Run OA

        sh /work/users/s/e/seyoun/CQTL_sQTL/scripts/namechange.sh {input.VCF} {input.OA} > {output.oa_matched}
        
        """

rule reheader:
    input:
        VCF_in=rules.create_vcf.output.vcf,
        matched_pbs=rules.processVCF.output.pbs_matched,
        matched_fnf=rules.processVCF.output.fnf_matched,
        matched_oa=rules.processVCF.output.oa_matched
    output:
        pbs_vcf="output/geno/pbs_geno/pbs.rename_wCHR.vcf.gz",
        fnf_vcf="output/geno/fnf_geno/fnf.rename_wCHR.vcf.gz",
        oa_vcf="output/geno/oa_geno/oa.rename_wCHR.vcf.gz"
    params:
        samtoolsVer=config['samtoolsVers']
        
    shell:
        """
        ml samtools/{params.samtoolsVer}

        bcftools reheader --sample {input.matched_pbs} {input.VCF_in} | bcftools annotate --rename-chrs output/geno/chrnm.txt -Oz -o {output.pbs_vcf}
        tabix -p vcf {output.pbs_vcf}
        
        bcftools reheader --sample {input.matched_fnf} {input.VCF_in} | bcftools annotate --rename-chrs output/geno/chrnm.txt -Oz -o {output.fnf_vcf}
        tabix -p vcf {output.fnf_vcf}
        bcftools reheader --sample {input.matched_oa} {input.VCF_in} | bcftools annotate --rename-chrs output/geno/chrnm.txt -Oz -o {output.oa_vcf}
        tabix -p vcf {output.oa_vcf}

        
        """

