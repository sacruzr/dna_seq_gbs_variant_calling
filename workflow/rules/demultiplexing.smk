rule get_barcodemap:
    output:
        barcodemap = "results/{plate}/reads/barcodemap.tsv"
    input:
        unpack(get_demultiplex_fastqs)
    run:
        with open(output.barcodemap, "w") as f:
            for n, row in samples[samples['plate'] == wildcards.plate].iterrows():
                linde_id = row['line_id']
                barcode = row['barcode']
                print(barcode, linde_id, sep="\t", file=f)
                
                
rule seqtk_subsample_pe:
    input:
        unpack(get_subsample_fastqs)
    output:
        f1="data/{plate}.1.fastq.gz",
        f2="data/{plate}.2.fastq.gz"
    params:
        n=500000,
        seed=12345
    threads:
        1
    wrapper:
        "v1.21.1/bio/seqtk/subsample/pe"

checkpoint demultiplex_pe:
    output: 
        outdir = temp(directory("results/{plate}/reads/demultiplexing")),
        logfile = "results/{plate}/logs/process_radtags.data.log"
    input: 
        unpack(get_demultiplex_fastqs),
        barcodemap = "results/{plate}/reads/barcodemap.tsv",
    params:
        extra = "-e {enzyme}".format(enzyme = config['demultiplexing']['re_enzime'])
    threads:
        resources['stacks']['threads']
    log:
        "results/{plate}/logs/demultiplexing.log"
    conda:
        "demultiplexing"
    shell:
        """
        mkdir {output.outdir}
        process_radtags -P \
            --threads {threads} \
            -1 {input.r1} \
            -2 {input.r2} \
            -b {input.barcodemap} \
            -o {output.outdir} \
            --retain-header \
            --inline_null \
            -c -r \
            {params.extra} 2> {log} &&
        mv {output.outdir}/process_radtags.data.log {output.logfile}
        """

rule concat_remaining_reads:
    input:
        get_remaining_reads
    output:
        temp("results/{plate}/reads/demultiplexing/remaining/{sample}.rem.fq.gz")
    shell:
        """
        zcat {input} | gzip -c > {output}
        """