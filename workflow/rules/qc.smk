rule fastqc:
    input:
        unpack(get_sample_fastqc),
    output:
        html="results/{plate}/qc/fastqc/{sample}.html",
        zip="results/{plate}/qc/fastqc/{sample}.zip",
    params: "--quiet"
    threads: resources['fastqc']['sample']
    wrapper:
        "v1.21.1/bio/fastqc"
        
rule fastqc_library:
    input:
        unpack(get_library_fastqc),
    output:
        html="results/{plate}/qc/fastqc_run/{sq_unit}.html",
        zip="results/{plate}/qc/fastqc_run/{sq_unit}.zip",
    params: "--quiet"
    log:
        "results/{plate}/logs/fastqc/{sq_unit}.log",
    threads: resources['fastqc']['library']
    wrapper:
        "v1.21.1/bio/fastqc"


rule samtools_stats:
    input:
        bam="results/{plate}/mapping/{sample}.sorted.bam"
    output:
        "results/{plate}/mapping/samtools-stats/{sample}.txt"
    wrapper:
        "v1.21.1/bio/samtools/stats"
        
rule qual_stats_read_pos:
    input:
        bam="results/{plate}/mapping/{sample}.sorted.bam",
        ref="resources/genome.fasta",
    output:
        "results/{plate}/mapping/readpos_stats/{sample}_readpos.stats"
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
        merged_file = "results/{plate}/stats/readpos/plateQC_stats.csv",
        plot_file = "results/{plate}/stats/readpos/plateQC_readpos.pdf",
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
        plate_qc = 'results/{plate}/qc/plateQC.json'
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
        
        with open(output.plate_qc, "w") as write_file:
            json.dump(plate_qc_data, write_file, indent=4)
        
rule multiqc:
    input: 
        get_multiqc_files
    output:
        "results/{plate}/qc/multiqc.html"
    params: 
        use_input_files_only=True
    log:
        "results/{plate}/logs/multiqc.log",
    wrapper:
        "v1.21.1/bio/multiqc"
