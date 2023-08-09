import pandas as pd
import os
import glob
from yaml import safe_load
import random
import string
import numpy as np
# Config file and sample sheet

configfile: "config/config.yaml"
    

# file containing sequencing unit info. flowcell, lane, index
# and corresponding plate number.
sequencing_units = pd.read_table(config["sequencing_runs"], sep='\t')
sequencing_units['plate'] = sequencing_units.plate.astype(str).str.pad(width = 2,side='left', fillchar='0')
sequencing_units.set_index('plate', drop=False, inplace=True)

references = pd.read_table(config["references_file"], sep='\t')
references.set_index('ref_name', drop=False, inplace=True)

samples = pd.read_table(config["barcode_map"], sep='\t')
samples = samples.merge(sequencing_units[['sequencing_unit_id', 'plate', 'fq1', 'fq2']], on='sequencing_unit_id')

#samples.set_index(["sample", "plate"], drop=False, inplace = True)
#samples.index = samples.index.set_levels(
#    [i.astype(str) for i in samples.index.levels]
#)  # enforce str in index





# resources params
with open(config["resources_config"], "r") as f:
    resources = safe_load(f)

def get_subsample_fastqs(wildcards):
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    if pd.isna(fastqs.fq2):
        print("multiple runs present in the same plate... UNEXPECTED BEHAVIOR")
    else:
        return {"f1": fastqs.fq1, "f2": fastqs.fq2}
   
    
def get_str_ref(wildcards):
    return references.loc[wildcards.ref]['STR_path']

def get_demultiplex_fastqs(wildcards):
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    
    if not pd.isna(fastqs.fq2):
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    else:
        return {"r1": fastqs.fq1}
        
def get_rawread_folder_name(wildcards):
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    raw_read_folder = fastqs.fq1.split('/')[-2]
    return raw_read_folder
    
def get_demultiplex_params(wildcards):
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    if not pd.isna(fastqs.fq2):
        return "-1 {R1} -2 {R2}".format(R1=fastqs.fq1, R2=fastqs.fq2)
    else:
        return "-f {R1}".format(R1=fastqs.fq1)
    
def get_remaining_reads(wildcards):
    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output.outdir
    
    return glob.glob(checkpoint_output + "/{sample}.rem.[1,2].fq.gz".format(sample=wildcards.sample))

def get_concat_remaining_reads(wildcards):
    return "results/{plate}/reads/demultiplexing/remaining/{sample}.rem.fq.gz".format(plate=wildcards.plate, sample=wildcards.sample)

def get_remaining_unpaired_trimmed_reads(wildcards):
    rem = "results/{plate}/reads/trimmed/paired/remaining/{sample}.fastq.gz".format(plate=wildcards.plate, sample=wildcards.sample)
    unpaired = glob.glob("results/{plate}/reads/trimmed/paired/{sample}.[1,2].unpaired.fastq.gz".format(plate=wildcards.plate, sample=wildcards.sample))
    unpaired.append(rem)
    return unpaired
def get_concat_remaining_unpaired_reads(wildcards):
    return "results/{plate}/reads/trimmed/paired/{sample}.rem.unpaired.fastq.gz".format(plate=wildcards.plate, sample=wildcards.sample)

def get_bams_to_merge(wildcards):
    se = "results/{plate}/mapping/{ref}/bwa/rem_unpaired/{sample}.sorted.bam".format(plate=wildcards.plate,
                                                                                     sample=wildcards.sample,
                                                                                    ref = wildcards.ref)
    pe = "results/{plate}/mapping/{ref}/bwa/mapping/{sample}.sorted.bam".format(plate=wildcards.plate,
                                                                                sample=wildcards.sample,
                                                                               ref = wildcards.ref)
    return [se, pe]

def get_bams_to_sort(wildcards):
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    if not pd.isna(fastqs.fq2):
        bam = "results/{plate}/mapping/{ref}/bwa/{sample}_merged.sorted.bam".format(**wildcards)
        return  bam
    else:
        bam = "results/{plate}/mapping/{ref}/bwa/mapping/{sample}.sorted.bam".format(**wildcards)
        return bam


def get_multiqc_files(wildcards):
    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output.outdir
    
    plate_df = sequencing_units.loc[wildcards.plate]
    samples_df = samples[samples['sequencing_unit_id'] == plate_df['sequencing_unit_id']]
    sample_names = samples_df['line_id'].tolist()
    multiqc_files = expand([
        "results/{plate}/mapping/{ref}/samtools-stats/{sample}.txt"], sample = sample_names, plate = wildcards.plate, ref = wildcards.ref)
    
    return multiqc_files

# function to generate raw reads name for fastqc_library execution
def get_sequencing_units_qc(wildcards):
    fastqs = list()
    
    for n, row in sequencing_units.iterrows():
        
        outs = expand([
            "results/{plate}/qc/fastqc_run/{sq_unit}-{group}.html",
            "results/{plate}/qc/fastqc_run/{sq_unit}-{group}.zip",
            
        ],plate = row.plate, sq_unit = row.sequencing_unit_id, group = ['R1', 'R2'])
        fastqs.extend(outs)
    return fastqs

def get_plate_efficiencies(wildcards):
    
    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output.outdir

    sample_list = glob.glob(checkpoint_output + "/*[!rem]*.fq.gz")
    sample_names = list(set([s.split('/')[-1].split('.')[0] for s in sample_list]))
    
    samtools_stats = expand([
        "results/{plate}/mapping/{ref}/samtools-stats/{sample}.txt",
    ], sample = sample_names, plate = wildcards.plate, ref=wildcards.ref)
    
    stacks_log = "results/{plate}/logs/process_radtags.data.log".format(plate = wildcards.plate)

    return {"samtools_stats": samtools_stats, "stacks_log": stacks_log}
        



def get_library_fastqc(wildcards):
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    if wildcards.group == 'R1':
        return fastqs.fq1
    elif wildcards.group == 'R2':
        return fastqs.fq2
    else:
        print("UNEXPECTED BEHAVIOR")

    

def get_sample_fastq(wildcards):
    
    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output.outdir
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    
    if not pd.isna(fastqs.fq2):
        sample_fastq = glob.glob(checkpoint_output + "/{sample}.[1,2].fq.gz".format(sample = wildcards.sample))
        return {"r1": sample_fastq[0], "r2": sample_fastq[1]}
    else:
        sample_fastq = glob.glob(checkpoint_output + "/{sample}.fq.gz".format(sample = wildcards.sample))
        return sample_fastq
        
def get_trimmed_reads(wildcards):
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    if not pd.isna(fastqs.fq2):
        return expand(
            "results/{plate}/reads/trimmed/paired/{sample}.{group}.fastq.gz",
            group = [1,2],
            **wildcards
        )
    else:
        return expand(
            "results/{plate}/reads/trimmed/single/{sample}.fastq.gz",
            **wildcards
        )
    

def get_readpos_files(wildcards):
    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output.outdir

    plate_df = sequencing_units.loc[wildcards.plate]
    samples_df = samples[samples['sequencing_unit_id'] == plate_df['sequencing_unit_id']]
    sample_names = samples_df['line_id'].tolist()
    readpos_files = expand([
        "results/{plate}/mapping/{ref}/readpos_stats/{sample}_readpos.stats",
    ], sample = sample_names, plate = wildcards.plate, ref = wildcards.ref)
    return readpos_files

wildcard_constraints:
    sq_unit = "|".join(sequencing_units['sequencing_unit_id'].unique()),
    plate="|".join(sequencing_units['plate'].unique()),
    sample="|".join(samples['line_id'].unique()),
    population="|".join(samples['population_name'].unique()),
    ref = "|".join(references['ref_name'].unique()),
    maf = "|".join(['0.05']),
    dp = "|".join([str(d) for d in range(1,20)]),
    qual = "|".join([str(q) for q in range(1,60)]),
    
def get_big_temp(wildcards):
    """Sets a temp dir for rules that need more temp space that is typical on some cluster environments. Defaults to system temp dir."""
    if config['bigtmp']:
        if config['bigtmp'].endswith("/"):
            return config['bigtmp'] + "".join(random.choices(string.ascii_uppercase, k=12)) + "/"
        else:
            return config['bigtmp'] + "/" + "".join(random.choices(string.ascii_uppercase, k=12)) + "/"
    else:
        return tempfile.gettempdir()
    
def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform="ILLUMINA",
    )

def get_mapped_reads(wildcards):
    return expand([
        "results/{plate}/mapping/{ref}/bwa/mapping/{sample}.sorted.bam",
        ],**wildcards)

def first_variant_calling():
    vcfs = list()
    for ref in list(references.ref_name.unique()):
        for plate, row in samples.iterrows():
            i_vcf = "results/{plate}/mapping/{ref}/first_variant_calling/{sample}_bwa_NGSEP.vcf.gz".format(plate = row.plate,
                                                                                     ref = ref,
                                                                                     sample = row.line_id)
            vcfs.append(i_vcf)
    return vcfs

def snps_plate():
    refs = references.ref_name.unique()
    plates = sequencing_units['plate'].tolist()
    quality = config['plate_snps']['Q']
    dp = config['plate_snps']['dp']
    maf = config['plate_snps']['maf']
    summaries = expand('results/{plate}/mapping/{ref}/vcf/{plate}_merged_Q{qual}_Dp{dp}_MAF{maf}_sorted_summary.txt',
        plate = plates,
        ref = refs,
        qual = quality,
        dp = dp,
        maf = maf)
    return summaries

def get_seq_ref_list(wildcards):
    ref_tab = pd.read_csv('resources/genome.fasta.fai', sep='\t', header=None)
    return ref_tab[0].tolist()

def get_sample_vcfs_by_population_merge_variants(wildcards):
    pop_df = samples[samples['population_name'].str.contains(wildcards.population)]
    vcfs = list()
    for n, row in pop_df.iterrows():
        i_vcf = 'results/{plate}/mapping/{ref}/{sample}_bwa_NGSEP.vcf.gz'.format(plate=row.plate, 
                                                                                 sample=row.line_id,
                                                                                 ref=wildcards.ref)
        vcfs.append(i_vcf)
    return vcfs

def get_sample_vcfs_by_plate_merge_variants(wildcards):
    plate_df = samples[samples['plate'] == wildcards.plate]
    vcfs = list()
    for n, row in plate_df.iterrows():
        i_vcf = 'results/{plate}/mapping/{ref}/first_variant_calling/{sample}_bwa_NGSEP.vcf.gz'.format(plate=row.plate, 
                                                                                 sample=row.line_id,
                                                                                 ref=wildcards.ref)
        vcfs.append(i_vcf)
    return vcfs

def get_sample_vcfs_by_plate_merge_vcfs(wildcards):
    plate_df = samples[samples['plate'] == wildcards.plate]
    vcfs = list()
    for n, row in plate_df.iterrows():
        i_vcf = "results/{plate}/mapping/{ref}/vcf/second_variant_call_plate/{sample}_bwa_NGSEP.vcf.gz".format(plate=wildcards.plate,
                                                                                                      sample=row.line_id,
                                                                                                     ref = wildcards.ref)
        vcfs.append(i_vcf)
    return vcfs
def get_sample_vcfs_by_population_merge_vcfs(wildcards):
    pop_df = samples[samples['population_name'].str.contains(wildcards.population)]
    vcfs = list()
    for n, row in pop_df.iterrows():
        i_vcf = "results/population/{population}/{ref}/vcfs/individual/{sample}_bwa_NGSEP.vcf".format(population=wildcards.population,
                                                                                                      sample=row.line_id,
                                                                                                     ref = wildcards.ref)
        vcfs.append(i_vcf)
    return vcfs


def get_bam_sample_population(wildcards):
    pop_df = samples[samples['population_name'].str.contains(wildcards.population)]
    sample_df = pop_df[pop_df['line_id'] == wildcards.sample]
    bam = "results/{plate}/mapping/{ref}/{sample}.sorted.bam".format(plate = list(sample_df.plate.unique())[0],
                                                                     sample = wildcards.sample,
                                                                    ref = wildcards.ref)
    return bam