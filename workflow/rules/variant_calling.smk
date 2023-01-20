rule single_sample_variant_detector:
    input:
        bam="results/{plate}/mapping/{sample}.sorted.bam",
        ref="resources/genome.fasta",
    output:
        "results/{plate}/mapping/{sample}_bwa_NGSEP.vcf.gz"
    params:
        params = " ".join(["-{param} {value}".format(param=i_param, value = config['NGSEP']['SingleSampleVariantsDetector'][i_param]) for i_param in config['NGSEP']['SingleSampleVariantsDetector'].keys()]),
        mem = "-Xmx3g"
    log:
        "results/{plate}/mapping/{sample}_bwa_NGSEP.log"
    shell:
        """
        java {params.mem} -jar {config[NGSEP][path]} \
            SingleSampleVariantsDetector {params.params} \
            -r {input.ref} \
            -i {input.bam} \
            -o {output} 2> {log} &&
            gzip -c {output}.vcf > {output} &&
            rm {output}.vcf
        """
