rule first_filtering:
    input:
        'results/{plate}/mapping/{ref}/vcf/{plate}_merged.vcf.gz'
    output:
        'results/{plate}/mapping/{ref}/vcf/{plate}_merged_Q{qual}.vcf.gz'
    log:
        'results/{plate}/mapping/{ref}/vcf/{plate}_merged_Q{qual}.log'
    params:
        mem = "-Xmx40g"
    shell:
        """
        java {params.mem} -jar {config[NGSEP][path]} \
            VCFFilter -i {input} -q {wildcards.qual} 2> {log}| \
            {config[bgzip]} > {output}
        """
rule sort_plate_vcfs:
    input:
        vcf = "results/{plate}/mapping/{ref}/vcf/{plate}_merged_Q{qual}.vcf.gz"
    output:
        vcf = "results/{plate}/mapping/{ref}/vcf/{plate}_merged_Q{qual}_sorted.vcf.gz",
        index = "results/{plate}/mapping/{ref}/vcf/{plate}_merged_Q{qual}_sorted.vcf.gz.tbi"
    resources:
        tmpdir = get_big_temp
    shell:
        """
        if [[ ! -f "{input.vcf}.tbi" ]]
        then
            {config[tabix]} -p vcf {input.vcf}
        fi
        {config[bcftools]} sort -T {resources.tmpdir} -Oz -o {output.vcf} {input.vcf} && \
        {config[tabix]} -p vcf {output.vcf}
        """
        
rule snps_plate_vcfs:
    input:
        vcf = "results/{plate}/mapping/{ref}/vcf/{plate}_merged_Q{qual}_sorted.vcf.gz",
    output:
        snps_vcf = "results/{plate}/mapping/{ref}/vcf/{plate}_merged_Q{qual}_Dp{dp}_MAF{maf}_sorted.vcf.gz",
    
    resources:
        tmpdir = get_big_temp
    threads: 30
    shell:
        """
        {config[bcftools]} view -v snps -m 2 -M 2 {input.vcf} | \
        grep -v scaffold | \
        {config[bcftools]} +fill-tags | \
        {config[bcftools]} filter -Ou -S . -i '(FORMAT/GQ)>={wildcards.qual} & (FORMAT/DP)>={wildcards.dp} & MAF >={wildcards.maf}' | \
        {config[bgzip]} > {output.snps_vcf}
        """