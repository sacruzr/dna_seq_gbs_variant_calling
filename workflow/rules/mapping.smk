rule map_reads:
    input:
        reads = get_trimmed_reads,
        idx = rules.bwa_index.output,
    output:
        temp("results/{plate}/mapping/bwa/paired/{sample}.sorted.bam")
    log:
        "results/{plate}/logs/bwa/{sample}.log"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate",
    resources:
        tmpdir = get_big_temp
    threads: resources['bwa_mem']['threads']
    wrapper:
        "v1.21.1/bio/bwa/mem"              

rule map_reads_se:
    input:
        reads = get_concat_remaining_unpaired_reads,
        idx = rules.bwa_index.output,
    output:
        temp("results/{plate}/mapping/bwa/rem_unpaired/{sample}.sorted.bam")
    log:
        "results/{plate}/logs/bwa/rem_unpaired/{sample}.log"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate",
    resources:
        tmpdir = get_big_temp
    threads: resources['bwa_mem']['threads']
    wrapper:
        "v1.21.1/bio/bwa/mem"
        
rule merge_bams:
    input:
        get_bams_to_merge
    output:
        temp("results/{plate}/mapping/bwa/{sample}_merged.sorted.bam")
    threads: 
        resources['samtools']['merge']['threads']
    wrapper:
        "v1.21.1/bio/samtools/merge"

rule mark_duplicates:
    input:
        bams="results/{plate}/mapping/bwa/{sample}_merged.sorted.bam",
    output:
        bam=temp("results/{plate}/mapping/{sample}.bam"),
        metrics="results/{plate}/mapping/dedup/{sample}.metrics.txt",
    log:
        "results/{plate}/logs/picard/dedup/{sample}.log",
    params:
        config["picard"]["MarkDuplicates"],
    resources:
        mem_mb=resources['picard_markduplicateds']['mem'],
        tmpdir = get_big_temp
    wrapper:
        "v1.21.1/bio/picard/markduplicates"

rule sort_bam:
    input:
        "results/{plate}/mapping/{sample}.bam",
    output:
        "results/{plate}/mapping/{sample}.sorted.bam",
    params:
        extra=resources['samtools']['sort']['mem']
    resources:
        tmpdir = get_big_temp
    threads:
        resources['samtools']['sort']['threads']
    wrapper:
        "v1.21.1/bio/samtools/sort"
        
rule index_bam:
    input:
        "results/{plate}/mapping/{sample}.sorted.bam",
    output:
        "results/{plate}/mapping/{sample}.sorted.bam.bai",
    threads:
        resources['samtools']['index']['threads']
    wrapper:
        "v1.21.1/bio/samtools/index"
