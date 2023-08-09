rule bam2gvcf:
    input:
        ref="resources/genome.fasta",
        indexes=expand("resources/genome.fasta.{ext}", ext=["sa", "pac", "bwt", "ann", "amb", "fai"]),
        dictf="resources/genome.dict",
        bam="results/{plate}/mapping/dedup/{sample}.sorted.bam",
        bai="results/{plate}/mapping/dedup/{sample}.bam.bai",
        l="resources/intervals/gvcf_intervals/{l}-scattered.interval_list"
    output:
        gvcf = "results/{plate}/interval_gvcfs/{sample}/{l}.raw.g.vcf.gz",
        gvcf_idx = "results/{plate}/interval_gvcfs/{sample}/{l}.raw.g.vcf.gz.tbi"
    resources:
        #!The -Xmx value the tool is run with should be less than the total amount of physical memory available by at least a few GB
        # subtract that memory here
        mem_mb = lambda wildcards, attempt: attempt * resources['bam2gvcf']['mem'],   # this is the overall memory requested
        reduced = lambda wildcards, attempt: attempt * (resources['bam2gvcf']['mem']- 3000)  # this is the maximum amount given to java
    log:
        "results/{plate}/logs/gatk_hc/{sample}/{l}.txt"
    params:
        minPrun = config['minP'],
        minDang = config['minD'],
    conda:
        "../envs/bam2vcf.yaml"
    shell:
        "gatk HaplotypeCaller "
        "--java-options \"-Xmx{resources.reduced}m\" "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output.gvcf} "
        "-L {input.l} "
        "--emit-ref-confidence GVCF --min-pruning {params.minPrun} --min-dangling-branch-length {params.minDang} &> {log}"
    
rule concat_gvcfs:
    input:
        unpack(get_interval_gvcfs)
    output:
        gvcf = "results/{plate}/gvcfs/{sample}.g.vcf.gz",
        tbi = "results/{plate}/gvcfs/{sample}.g.vcf.gz.tbi"
    conda:
        "demultiplexing"
    resources:
        tmpdir = get_big_temp
    shell:
        """
        bcftools concat -D -a -Ou {input} | bcftools sort -T {resources.tmpdir} -Oz -o {output.gvcf} -
        tabix -p vcf {output.gvcf}
        """

rule create_db_mapfile:
    input:
        get_input_mapfile
    output:
        db_mapfile="results/{plate}/genomics_db_import/DB_mapfile.txt"
    run:
        with open(output.db_mapfile, "w") as f:
            for file_path in input:
                sample_name = os.path.basename(file_path).replace(".g.vcf.gz", "")
                print(sample_name, file_path, sep="\t", file=f)

rule prepare_db_intervals:
    """GenomicsDBImport needs list of intervals to operate on so this rule writes that file"""
    input:
        fai = "resources/genome.fasta.fai",
    output:
        intervals = "results/{plate}/genomics_db_import/db_intervals.list"
    run:
        with open(output.intervals, "w") as out:
            with open(input.fai, "r") as f:
                for line in f:
                    line = line.strip().split()
                    chrom, end = line[0], line[1]
                    print(f"{chrom}:1-{end}", file=out)

rule gvcf2DB:
    input:
        unpack(get_gvcfs_db),
        db_mapfile = "results/{plate}/genomics_db_import/DB_mapfile.txt",
        intervals = "results/{plate}/genomics_db_import/db_intervals.list"
    output:
        db = directory("results/{plate}/genomics_db_import/DB"),
        tar = "results/{plate}/genomics_db_import/DB.tar"
    log:
        "results/{plate}/logs/gatk_db_import.txt"
    conda:
        "../envs/bam2vcf.yaml"
    resources:
        #!The -Xmx value the tool is run with should be less than the total amount of physical memory available by at least a few GB
        # subtract that memory here
        mem_mb = lambda wildcards, attempt: attempt * resources['gvcf2DB']['mem'],   # this is the overall memory requested
        reduced = lambda wildcards, attempt: attempt * (resources['gvcf2DB']['mem']- 3000)  # this is the maximum amount given to java
    shell:
        # NOTE: reader-threads > 1 useless if you specify multiple intervals
        # a forum suggested TILEDB_DISABLE_FILE_LOCKING=1 to remedy sluggish performance
        """
        export TILEDB_DISABLE_FILE_LOCKING=1
        gatk GenomicsDBImport \
            --java-options '-Xmx{resources.reduced}m -Xms{resources.reduced}m' \
            --genomicsdb-shared-posixfs-optimizations true \
            --batch-size 25 \
            --genomicsdb-workspace-path {output.db} \
            -L {input.intervals} \
            --merge-input-intervals \
            --sample-name-map {input.db_mapfile} &> {log}

        tar -cf {output.tar} {output.db}
        """

rule DB2vcf:
    input:
        db = "results/{plate}/genomics_db_import/DB.tar",
        ref = "resources/genome.fasta",
    output:
        vcf = "results/{plate}/vcfs/raw.vcf.gz",
        vcfidx = "results/{plate}/vcfs/raw.vcf.gz.tbi"
    params:
        het = config['het_prior'],
        db = lambda wildcards, input: input.db[:-4]
    resources:
         #!The -Xmx value the tool is run with should be less than the total amount of physical memory available by at least a few GB
        # subtract that memory here
        mem_mb = lambda wildcards, attempt: attempt * resources['DB2vcf']['mem'],   # this is the overall memory requested
        reduced = lambda wildcards, attempt: attempt * (resources['DB2vcf']['mem']- 3000),  # this is the maximum amount given to java
        tmpdir = get_big_temp
    log:
        "result/{plate}/logs/gatk_genotype_gvcfs.txt"
    conda:
        "../envs/bam2vcf.yaml"
    shell:
        """
        tar -xf {input.db}
        gatk GenotypeGVCFs \
            --java-options '-Xmx{resources.reduced}m -Xms{resources.reduced}m' \
            -R {input.ref} \
            --heterozygosity {params.het} \
            --genomicsdb-shared-posixfs-optimizations true \
            -V gendb://{params.db} \
            -O {output.vcf} \
            --tmp-dir {resources.tmpdir} &> {log}
        """
