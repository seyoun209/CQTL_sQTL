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
