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
                
                

checkpoint demultiplex:
    output: 
        outdir = directory("results/{plate}/reads/demultiplexing/stacks"),
        logfile = "results/{plate}/logs/process_radtags.data.log"
    input: 
        unpack(get_demultiplex_fastqs),
        barcodemap = "results/{plate}/reads/barcodemap.tsv",
    params:
        extra = "-e {enzyme}".format(enzyme = config['demultiplexing']['re_enzime']),
        raw_reads_folder = lambda wildcards: get_rawread_folder_name(wildcards),
        files = lambda wildcards: get_demultiplex_params(wildcards)
    threads:
        resources['stacks']['threads']
    log:
        "results/{plate}/logs/demultiplexing.log"
    conda:
        "../envs/demultiplexing.yaml"
    shell:
        """
        mkdir {output.outdir}
        process_radtags \
            --threads {threads} \
            {params.files} \
            -b {input.barcodemap} \
            -o {output.outdir} \
            --retain-header \
            --inline_null \
            -r \
            {params.extra} 2> {log} &&
        mv {output.outdir}/process_radtags.{params.raw_reads_folder}.log {output.logfile}
        """


        
rule concat_remaining_reads:
    input:
        get_remaining_reads
    output:
        temp("results/{plate}/reads/demultiplexing/remaining/{sample}.rem.fq.gz")
    shell:
        """
        zcat {input} | {config[bgzip]} -c > {output}
        """