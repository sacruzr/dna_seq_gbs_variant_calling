rule genome_faidx:
    input:
        path = "resources/{ref}.fasta"
    output:
        "resources/{ref}.fasta.fai"
    cache: True
    wrapper:
        "v1.21.1/bio/samtools/faidx"

rule bwa_index:
    input:
        "resources/{ref}.fasta",
    output:
         idx=multiext("resources/{ref}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/{ref}_bwa_index.log"
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.21.1/bio/bwa/index"
        
rule list_ref_seqs:
    input:
        ref_index = "resources/{ref}.fasta.fai"
    output:
        out_list = 'resources/{ref}.seqlist'
    run:
        ref_tab = pd.read_csv(input.ref_index, sep='\t', header=None)
        ref_tab[0].to_csv(output.out_list, sep='\t', index =False)