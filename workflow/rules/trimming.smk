rule trim_reads_pe:
    input:
        unpack(get_sample_fastq)
    output:
        r1="results/{plate}/reads/trimmed/paired/{sample}.1.fastq.gz",
        r2="results/{plate}/reads/trimmed/paired/{sample}.2.fastq.gz",
        r1_unpaired=temp("results/{plate}/reads/trimmed/paired/{sample}.1.unpaired.fastq.gz"),
        r2_unpaired=temp("results/{plate}/reads/trimmed/paired/{sample}.2.unpaired.fastq.gz"),
        trimlog=temp("results/{plate}/reads/trimmed/paired/{sample}.trimlog.txt")
    params:
        **config["trimmomatic"]["pe"],
        extra=lambda w, output: "-trimlog {output}".format(output = output.trimlog)
    log:
        "results/{plate}/logs/trimmomatic/{sample}.log"
    threads:
        resources['trimmomatic']['threads']
    resources:
        mem_mb=resources['trimmomatic']['mem']
    wrapper:
        "v1.21.1/bio/trimmomatic/pe"
        
rule trim_reads_se:
    input:
        get_sample_fastq,
    output:
        "results/{plate}/reads/trimmed/single/{sample}.fastq.gz"
    params:
        **config["trimmomatic"]["pe"]
    log:
        "results/{plate}/logs/trimmomatic/{sample}.log"
    threads:
        resources['trimmomatic']['threads']
    resources:
        mem_mb=resources['trimmomatic']['mem']
    wrapper:
        "v1.21.1/bio/trimmomatic/se"
        
rule trim_remaining_reads:
    input:
        get_concat_remaining_reads
    output:
        temp("results/{plate}/reads/trimmed/paired/remaining/{sample}.fastq.gz")
    params:
        **config["trimmomatic"]["pe"]
    log:
        "results/{plate}/logs/trimmomatic/remaining/{sample}.log"
    threads:
        resources['trimmomatic']['threads']
    resources:
        mem_mb=resources['trimmomatic']['mem']
    wrapper:
        "v1.21.1/bio/trimmomatic/se"
        
rule concat_remaning_unpaired_trimmed_reads:
    input:
        get_remaining_unpaired_trimmed_reads
    output:
        "results/{plate}/reads/trimmed/paired/{sample}.rem.unpaired.fastq.gz"
        
    shell:
        """
        zcat {input} | {config[bgzip]} -c > {output}
        """