rule multiqc:
    input:
        fastqc=[expand("output/QC/{sampleName}_{read}_fastqc.{ext}", sampleName=key, read=['R1', 'R2'], ext=['zip']) for key in read1],
        trimming=[expand('output/trim/{sampleName}_{ext}.fastq.gz_trimming_report.txt', sampleName=key, ext=['R1','R2']) for key in read1],
        alignment=expand("output/align/{sampleName}_Log.final.out", sampleName=read1.keys()),
        quantification=expand("output/quant/{sampleName}/qunat.sf", sampleName=read1.keys()),
        verifybamid=[expand('output/QC/{sampleName}_verifybamid.{ext}', sampleName=key, ext=['bestSM','selfSM']) for key in read1]
    output:
          multiqc=directory("output/multiqc")
    params:
        version = config['multiqcVers']
    log:
        out = 'output/logs/multiqc.out',
        err = 'output/logs/multiqc.err'
    shell:
        """
        ml multiqc/{params.version}
        mkdir -p output/multiqc
        multiqc -f -o {output.multiqc} {input.fastqc} {input.trimming} {input.alignment} {input.quantification} {input.verifybamid} 1> {log.out} 2> {log.err}

        """

