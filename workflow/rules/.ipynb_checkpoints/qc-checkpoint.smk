rule fastqc:
    input:
        unpack(get_sample_fastq),
    output:
        html="results/{plate}/qc/fastqc/{sample}.html",
        zip="results/{plate}/qc/fastqc/{sample}.zip",
    params: "--quiet"
    threads: resources['fastqc']['sample']
    wrapper:
        "v1.21.1/bio/fastqc"
        
rule fastqc_library:
    input:
        get_library_fastqc
    output:
        html="results/{plate}/qc/fastqc_run/{sq_unit}-{group}.html",
        zip="results/{plate}/qc/fastqc_run/{sq_unit}-{group}.zip",
    params: "--quiet"
    log:
        "results/{plate}/logs/fastqc/{sq_unit}_{group}.log",
    threads: resources['fastqc']['library']
    wrapper:
        "v1.21.1/bio/fastqc"


rule samtools_stats:
    input:
        bam="results/{plate}/mapping/{ref}/{sample}.sorted.bam"
    output:
        "results/{plate}/mapping/{ref}/samtools-stats/{sample}.txt"
    wrapper:
        "v1.21.1/bio/samtools/stats"
        
rule qual_stats_read_pos:
    input:
        bam="results/{plate}/mapping/{ref}/{sample}.sorted.bam",
        ref="resources/{ref}.fasta",
    output:
        "results/{plate}/mapping/{ref}/readpos_stats/{sample}_readpos.stats"
    params:
        mem = "-Xmx3g"
    shell:
        """
        java {params.mem} -jar {config[NGSEP][path]} BasePairQualStats \
            -o {output} -r {input.ref} {input.bam}
        """
    
rule plate_readpos_qc:
    input:
        readpos_stats = get_readpos_files
    output:
        merged_file = "results/{plate}/mapping/{ref}/stats/readpos/plateQC_stats.csv",
        plot_file = "results/{plate}/mapping/{ref}/stats/readpos/plateQC_readpos.pdf",
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import numpy as np
        
        formated_outs = list()
        for readpos_out_path in input.readpos_stats:
            # magic -14 refers to file extension "_readpos.stats"
            sample = readpos_out_path.split('/')[-1][:-14]
            # ignore last 3 lines
            readpos_out = pd.read_csv(readpos_out_path, skipfooter=3, sep='\t', header=None, engine='python')
            if readpos_out.shape[0] > 1:
                readpos_out['sample'] = sample
                readpos_out['mismatch_mm'] = (readpos_out[1]/readpos_out[3])*100
                readpos_out['mismatch_um'] = (readpos_out[2]/readpos_out[4])*100
                formated_outs.append(readpos_out)
            else:
                print("check sample {sample} seems to be empty the readpos stats file".format(sample = sample))
        merged_readpos = pd.concat(formated_outs, ignore_index = True)
        means = merged_readpos.groupby(0, as_index = False).mean(numeric_only=True)
        fig = plt.figure(figsize=(15,8))
        ax = fig.add_subplot(111)
        ax.plot(means[0], means['mismatch_mm'], label='Multialignments')
        ax.plot(means[0], means['mismatch_um'], label='Unique alignments')
        ax.set_xlabel('Read Position')
        ax.set_ylabel('Missmatch Percentage (/%)')
        xticks = list(np.arange(0, 151,2))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks, rotation=45, size=8)
        plt.grid(axis='y', color='0.9')
        plt.grid(axis='x', color='0.8')
        plt.legend()
        plt.savefig(output.plot_file)     
        merged_readpos.to_csv(output.merged_file, header = None, index = False)
        
rule plate_efficiencies_stats:
    input:
        unpack(get_plate_efficiencies)
    output:
        plate_qc = 'results/{plate}/mapping/{ref}/stats/plateQC.json'
    run:
        import os
        import json
        
        plate_qc_data = dict()
        # get total counts of stacks output
        data=dict()
        flag=False
        with open(input.stacks_log,'r') as f:
            for line in f:
                if line.startswith('BEGIN total_raw_read_counts'):
                    flag=True
                elif line.strip().endswith('END total_raw_read_counts'):
                    flag=False
                elif flag:
                    i_data = line[:-1].split('\t')
                    data[i_data[0]] = i_data[1:]
        
        plate_qc_data['total_reads'] = int(data['Total Sequences'][0])
        plate_qc_data['reads_not_found'] = plate_qc_data['total_reads'] - int(data['Retained Reads'][0])
        # aggregate mapping stats of all samples
        
        plate_qc_data['Gbp_non_adapters'] = 0
        plate_qc_data['Gbp_mapped'] = 0
        
        flag=False
        for i_file in input.samtools_stats:
            samtools_data = dict()
            with open(i_file,'r') as f:
                    for line in f:
                        if line.startswith('SN'):
                            flag=True
                            i_data = line[3:].split(':')
                            samtools_data[i_data[0]] = i_data[1]
                        elif not line.startswith('SN'):
                            flag=False
                            
            plate_qc_data['Gbp_non_adapters'] += int(samtools_data["total length"].split("\t")[1])
            plate_qc_data['Gbp_mapped'] += int(samtools_data["bases mapped (cigar)"].split("\t")[1])
        plate_qc_data['Gbp_non_adapters'] = plate_qc_data['Gbp_non_adapters']/1000000000
        plate_qc_data['Gbp_mapped'] = plate_qc_data['Gbp_mapped']/1000000000
        with open(output.plate_qc, "w") as write_file:
            json.dump(plate_qc_data, write_file, indent=4)
        
rule multiqc:
    input: 
        get_multiqc_files
    output:
        "results/{plate}/mapping/{ref}/stats/multiqc.html"
    params: 
        use_input_files_only=True
    log:
        "results/{plate}/logs/{ref}_multiqc.log",
    wrapper:
        "v1.21.1/bio/multiqc"
rule plate_mapping_stats:
    input:
        mapping_stats = 'results/{plate}/mapping/{ref}/stats/multiqc_data/multiqc_general_stats.txt'
    output:
        plate_stats_plot = 'results/{plate}/mapping/{ref}/stats/plate_stats_plot.pdf',
    run:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
        
        barcodes = pd.read_csv(config['barcodes'])
        
        mapping_stats = pd.read_csv(input.mapping_stats, sep='\t')
        mapping_stats.rename(columns = {'Sample': 'line_id'}, inplace = True)
        mapping_stats['log_total_reads'] = np.log(mapping_stats['Samtools_mqc-generalstats-samtools-raw_total_sequences']+1)
        mapping_responses = ['log_total_reads',
                            'Samtools_mqc-generalstats-samtools-reads_mapped_percent',
                            'Samtools_mqc-generalstats-samtools-error_rate']
        
        plate_data = samples.merge(barcodes, on = 'barcode')
        plate_data = plate_data.merge(mapping_stats, on = 'line_id' )
        
        
        fig, axs = plt.subplots(3, 1, figsize=(8, 10))
        for response, ax in zip(mapping_responses, axs):
            average = pd.pivot_table(data= plate_data, columns='well_r', index = 'well_c', values = response)
            sns.heatmap(average, cmap='viridis', ax = ax)
            ax.set_title(response)
        plt.suptitle(wildcards.plate + " ref:" + wildcards.ref)
        plt.tight_layout()
        plt.savefig(output.plate_stats_plot)
        
rule vcf_summary:
    input:
        snps_vcf = "results/{plate}/mapping/{ref}/vcf/{plate}_merged_Q{qual}_Dp{dp}_MAF{maf}_sorted.vcf.gz",
    output:
        summary = "results/{plate}/mapping/{ref}/vcf/{plate}_merged_Q{qual}_Dp{dp}_MAF{maf}_sorted_summary.txt",
    shell:
        """
        {config[bcftools]} stats -s - -d 0,100,2 {input.snps_vcf} > {output.summary}
        """

rule plate_snp_stats:
    input:
        snp_stats = "results/{plate}/mapping/{ref}/vcf/{plate}_merged_Q{qual}_Dp{dp}_MAF{maf}_sorted_summary.txt",
    output:
        plate_stats_plot = 'results/{plate}/mapping/{ref}/stats/plate_snp_Q{qual}_Dp{dp}_MAF{maf}_stats_plot.pdf',
    run:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
        
        barcodes = pd.read_csv(config['barcodes'])
        
        f = open(input.snp_stats, "r")
        lines = f.readlines()

        general_lines = [line.replace('\n', '').split('\t') for line in lines if re.match('^SN', line) != None]
        n_snps = int(general_lines[3][-1])

        taxa_lines = [line.replace('\n', '').split('\t') for line in lines if re.match('^PSC', line) != None]
        taxa_df = pd.DataFrame(taxa_lines)
        taxa_df.columns = ['PSC','id','line_id','nRefHom','nNonRefHom','nHets','nTransitions','nTransversions','nIndels','average_depth','nSingletons','nHapRef','nHapAlt','nMissing']
        taxa_df['total_variants'] = taxa_df['nRefHom'].astype(int) + taxa_df['nNonRefHom'].astype(int) + taxa_df['nHets'].astype(int)
        taxa_df['HO'] = taxa_df['nHets'].astype(int)/taxa_df['total_variants']
        taxa_df['TvTs'] =taxa_df['nTransitions'].astype(int) / taxa_df['nTransversions'].astype(int) 
        taxa_df['average_depth'] = taxa_df['average_depth'].astype(float)
        taxa_df['log_total_variants'] = np.log(taxa_df['total_variants']+1)
        plate_data = samples.merge(barcodes, on = 'barcode')
        plate_data = plate_data.merge(taxa_df, on = 'line_id' )
        
        
        snp_responses = ['log_total_variants', 'average_depth', 'HO', 'TvTs']
        
        fig, axs = plt.subplots(4, 1, figsize=(8, 10))
        for response, ax in zip(snp_responses, axs):
            print(response)
            average = pd.pivot_table(data= plate_data, columns='well_r', index = 'well_c', values = response)
            sns.heatmap(average, cmap='viridis', ax = ax)
            ax.set_title(response)
        plt.suptitle(wildcards.plate + " ref:" + wildcards.ref)
        plt.tight_layout()
        plt.savefig(output.plate_stats_plot)