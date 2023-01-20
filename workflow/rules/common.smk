import pandas as pd
import os
import glob
from yaml import safe_load
import random
import string

# Config file and sample sheet

configfile: "config/config.yaml"
    

# file containing sequencing unit info. flowcell, lane, index
# and corresponding plate number.
sequencing_units = pd.read_table(config["sequencing_runs"], sep='\t')
sequencing_units['plate'] = sequencing_units.plate.astype(str).str.pad(width = 2,side='left', fillchar='0')
sequencing_units.set_index('plate', drop=False, inplace=True)

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
    fastqs = samples[samples['plate'] == wildcards.plate][["fq1", "fq2"]].drop_duplicates(ignore_index=True).dropna()
    if fastqs.shape[0] > 1:
        print("multiple runs present in the same plate... UNEXPECTED BEHAVIOR")
    else:
        fastqs = fastqs.loc[0]
    if len(fastqs) == 2:
        return {"f1": fastqs.fq1, "f2": fastqs.fq2}
    return fastqs.fq1
    
    
#production function
#def get_demultiplex_fastqs(wildcards):
#    fastqs = samples[samples['plate'] == wildcards.plate][["fq1", "fq2"]].drop_duplicates(ignore_index=True).dropna()
#    if fastqs.shape[0] > 1:
#        print("multiple runs present in the same plate... UNEXPECTED BEHAVIOR")
#    else:
#        fastqs = fastqs.loc[0]
#    if len(fastqs) == 2:
#        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
#    return fastqs.fq1


#test function
def get_demultiplex_fastqs(wildcards):
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]

    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    else:
        print("UNEXPECTED BEHAVIOR")

def get_remaining_reads(wildcards):
    return glob.glob("results/{plate}/reads/demultiplexing/{sample}.rem.[1,2].fq.gz".format(plate=wildcards.plate, sample=wildcards.sample))

def get_concat_remaining_reads(wildcards):
    return "results/{plate}/reads/demultiplexing/remaining/{sample}.rem.fq.gz".format(plate=wildcards.plate, sample=wildcards.sample)

def get_remaining_unpaired_trimmed_reads(wildcards):
    rem = "results/{plate}/reads/trimmed/remaining/{sample}.fastq.gz".format(plate=wildcards.plate, sample=wildcards.sample)
    unpaired = glob.glob("results/{plate}/reads/trimmed/{sample}.[1,2].unpaired.fastq.gz".format(plate=wildcards.plate, sample=wildcards.sample))
    unpaired.append(rem)
    return unpaired
def get_concat_remaining_unpaired_reads(wildcards):
    return "results/{plate}/reads/trimmed/{sample}.rem.unpaired.fastq.gz".format(plate=wildcards.plate, sample=wildcards.sample)

def get_bams_to_merge(wildcards):
    se = "results/{plate}/mapping/bwa/rem_unpaired/{sample}.sorted.bam".format(plate=wildcards.plate, sample=wildcards.sample)
    pe = "results/{plate}/mapping/bwa/paired/{sample}.sorted.bam".format(plate=wildcards.plate, sample=wildcards.sample)
    return [se, pe]

def get_multiqc_files(wildcards):
    checkpoint_output = checkpoints.demultiplex_pe.get(**wildcards).output.outdir
    sample_list = glob.glob(checkpoint_output + "/*[!rem].1.fq.gz")
    sample_names = [s.split('/')[-1][:-8] for s in sample_list]
    multiqc_files = expand([
        "results/{plate}/mapping/samtools-stats/{sample}.txt",
        "results/{plate}/qc/fastqc/{sample}.zip",
        "results/{plate}/mapping/dedup/{sample}.metrics.txt"
    ], sample = sample_names, plate = wildcards.plate)
    
    return multiqc_files

# function to generate raw reads name for fastqc_library execution
def get_sequencing_units_qc():
    fastqs = list()
    for n, row in sequencing_units.iterrows():
        html_file= "results/{plate}/qc/fastqc_run/{sq_unit}.html".format(plate = row.plate, sq_unit = row.sequencing_unit_id)
        zip_file = "results/{plate}/qc/fastqc_run/{sq_unit}.zip".format(plate = row.plate, sq_unit = row.sequencing_unit_id)
        plate_qc_file = "results/{plate}/qc/plateQC.json".format(plate = row.plate)
        plate_qc_readpos_file = "results/{plate}/stats/readpos/plateQC_readpos.pdf".format(plate = row.plate)
        fastqs.extend([html_file, zip_file,plate_qc_file, plate_qc_readpos_file])
    return fastqs

def get_plate_efficiencies(wildcards):
    
    checkpoint_output = checkpoints.demultiplex_pe.get(**wildcards).output.outdir

    sample_list = glob.glob(checkpoint_output + "/*[!rem].1.fq.gz")
    sample_names = [s.split('/')[-1][:-8] for s in sample_list]
    
    samtools_stats = expand([
        "results/{plate}/mapping/samtools-stats/{sample}.txt",
    ], sample = sample_names, plate = wildcards.plate)
    
    stacks_log = "results/{plate}/logs/process_radtags.data.log".format(plate = wildcards.plate)

    return {"samtools_stats": samtools_stats, "stacks_log": stacks_log}
        



def get_library_fastqc(wildcards):
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    else:
        print("UNEXPECTED BEHAVIOR")

    

def get_sample_fastqc(wildcards):
    checkpoint_output = checkpoints.demultiplex_pe.get(**wildcards).output.outdir
    fastqs = glob.glob(checkpoint_output + "/{sample}.[1,2].fq.gz".format(sample = wildcards.sample))
    if len(fastqs) == 2:
        return {"r1": fastqs[0], "r2": fastqs[1]}
    else:
        print("UNEXPECTED BEHAVIOR")
        
def get_trimmed_reads(wildcards):
    return expand(
        "results/{plate}/reads/trimmed/{sample}.{group}.fastq.gz",
        group = [1,2],
        **wildcards
    )

def get_readpos_files(wildcards):
    checkpoint_output = checkpoints.demultiplex_pe.get(**wildcards).output.outdir
    sample_list = glob.glob(checkpoint_output + "/*[!rem].1.fq.gz")
    sample_names = [s.split('/')[-1][:-8] for s in sample_list]
    
    readpos_files = expand([
        "results/{plate}/mapping/readpos_stats/{sample}_readpos.stats",
    ], sample = sample_names, plate = wildcards.plate)
    return readpos_files

wildcard_constraints:
    sq_unit = "|".join(sequencing_units['sequencing_unit_id'].unique()),
    plate="|".join(sequencing_units['plate'].unique()),
    sample="|".join(samples['line_id'].unique()),
    
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
        "results/{plate}/mapping/bwa/paired/{sample}.sorted.bam",
        ],**wildcards)
