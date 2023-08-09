rule map_reads:
    input:
        reads = get_trimmed_reads,
        idx = rules.bwa_index.output,
    output:
        temp("results/{plate}/mapping/{ref}/bwa/mapping/{sample}.sorted.bam")
    log:
        "results/{plate}/logs/bwa/{ref}/{sample}.log"
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
        temp("results/{plate}/mapping/{ref}/bwa/rem_unpaired/{sample}.sorted.bam")
    log:
        "results/{plate}/logs/bwa/{ref}/rem_unpaired/{sample}.log"
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
        temp("results/{plate}/mapping/{ref}/bwa/{sample}_merged.sorted.bam")
    threads: 
        resources['samtools']['merge']['threads']
    wrapper:
        "v1.21.1/bio/samtools/merge"


rule sort_bam:
    input:
        get_bams_to_sort,
    output:
        "results/{plate}/mapping/{ref}/{sample}.sorted.bam",
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
        "results/{plate}/mapping/{ref}/{sample}.sorted.bam",
    output:
        "results/{plate}/mapping/{ref}/{sample}.sorted.bam.bai",
    threads:
        resources['samtools']['index']['threads']
    wrapper:
        "v1.21.1/bio/samtools/index"
