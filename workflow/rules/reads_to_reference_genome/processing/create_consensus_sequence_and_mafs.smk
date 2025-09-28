# Rule: Create consensus sequence using ANGSD
rule create_consensus_sequence:
    input:
        sorted_bam="{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.rescaled.bam",
        bam_index="{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.rescaled.bam.bai",
        reference_genome="{species}/raw/ref_genome/{ref_genome}.fa"
    output:
        consensus_sequence="{species}/processed/{ref_genome}/consensus/{individual}_{ref_genome}/{individual}_{ref_genome}_consensus.fa.gz"
    message: "Creating consensus sequence for {wildcards.individual} mapped to {wildcards.ref_genome} in species {wildcards.species}"
    log:
        "{species}/logs/{ref_genome}/consensus/{individual}/{individual}_{ref_genome}_consensus.log"
    conda:
        "../../../envs/angsd.yaml"
    shell:
        """
        angsd \
        -out {wildcards.species}/processed/{wildcards.ref_genome}/consensus/{wildcards.individual}_{wildcards.ref_genome}/{wildcards.individual}_{wildcards.ref_genome}_consensus  \
        -i {input.sorted_bam}  \
        -ref {input.reference_genome}  \
        -doFasta 2  \
        -doCounts 1  \
        -minMapQ 10  \
        -minQ 0  \
        -remove_bads 1  \
        -baq 1  \
        -C 50  \
        > {log} 2>&1
        """


# Rule: Create SNP and MAF files using ANGSD
rule create_snp:
    input:
        sorted_bam="{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.rescaled.bam",
        bam_index="{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.rescaled.bam.bai",
        reference_genome="{species}/raw/ref_genome/{ref_genome}.fa"
    output:
        mafs="{species}/processed/{ref_genome}/snp/{individual}_{ref_genome}/{individual}_{ref_genome}_snp.mafs.gz"
    message: "Creating SNP and MAF files for {wildcards.individual} mapped to {wildcards.ref_genome} in species {wildcards.species}"
    log:
        "{species}/logs/{ref_genome}/snp/{individual}/{individual}_{ref_genome}_snp.log"
    conda:
        "../../../envs/angsd.yaml"
    shell:
        """
        angsd \
        -out {wildcards.species}/processed/{wildcards.ref_genome}/snp/{wildcards.individual}_{wildcards.ref_genome}/{wildcards.individual}_{wildcards.ref_genome}_snp  \
        -i {input.sorted_bam}  \
        -ref {input.reference_genome}  \
        -doMaf 1  \
        -doMajorMinor 1  \
        -GL 2 \
        -SNP_pval 1e-6 \
        -minMapQ 10  \
        -minQ 0  \
        -remove_bads 1  \
        -baq 1  \
        -C 50  \
        > {log} 2>&1 
        """