rule get_ref:
    input:
        path = config['ref']['path']
    output:
        "resources/genome.fasta"
    shell:
        "zcat {input.path} > {output}"

rule genome_faidx:
    input:
        path = "resources/genome.fasta"
    output:
        "resources/genome.fasta.fai"
    cache: True
    wrapper:
        "v1.21.1/bio/samtools/faidx"

rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
         idx=multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index.log"
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.21.1/bio/bwa/index"
#
#rule genome_dict:
#    input:
#        "resources/genome.fasta"
#    output:
#        "resources/genome.dict"
#    conda:
#        "demultiplexing"
#    cache: True
#    shell:
#        "samtools dict {input} > {output}"
        #
#
#rule picard_intervals:
#    input:
#        ref="resources/genome.fasta",
#        fai="resources/genome.fasta.fai",
#        dictf="resources/genome.dict"
#    output:
#        intervals=temp("resources/intervals/picard_interval_list.list")
#    params:
#        minNmer = int(config['intervals']['minNmer'])
#    conda:
#        "../envs/bam2vcf.yaml"
#    log:
#        "logs/picard_intervals/log.txt"
#    shell:
        #"picard ScatterIntervalsByNs REFERENCE={input.ref} OUTPUT={output.intervals} MAX_TO_MERGE={params.minNmer} OUTPUT_TYPE=ACGT &> {log}"
#
#rule format_interval_list:
#    input:
#        intervals="resources/intervals/picard_interval_list.list"
#    output:
#        intervals="resources/intervals/master_interval_list.list"
#    run:
#        with open(output.intervals, "w") as out:
#            with open(input.intervals, "r") as inp:
#                for line in inp:
#                    if not line.startswith("@"):
#                        line = line.strip().split('\t')
#                        chrom, start, end = line[0], line[1], line[2]
#                        print(f"{chrom}:{start}-{end}", file=out)
#
#checkpoint create_gvcf_intervals:
#    input:
#        ref="resources/genome.fasta",
#        intervals="resources/intervals/master_interval_list.list"
#    output:
#        out_dir=directory("resources/intervals/gvcf_intervals")
#    params:
#        max_intervals=config["intervals"]["num_gvcf_intervals"]
#    log:
#        "logs/gvcf_intervals/log.txt"
#    conda:
#        "../envs/bam2vcf.yaml"
#    shell:
#        """
#        gatk SplitIntervals -L {input.intervals} \
#        -O {output} -R {input.ref} -scatter {params} \
#        -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
#        --interval-merging-rule OVERLAPPING_ONLY  &> {log}
#        """
