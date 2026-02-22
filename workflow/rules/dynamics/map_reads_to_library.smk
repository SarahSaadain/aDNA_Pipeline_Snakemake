
####################################################
# Snakemake rules
####################################################

# This rule maps sequencing reads of each individual to the combined SCG and TE library
rule map_reads_to_scg_feature_library:
    input:
        reads=["{species}/processed/reads/reads_merged/{individual}.fastq.gz"],
        # The index is produced by the bwa_index_library rule
        idx=multiext("{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"),
    output:
        temp("{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.sam"), # became temp to save space
    log:
        "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg_bwa.log",
    message: "Mapping reads of {wildcards.individual} to {wildcards.species} SCG and Feature library"
    params:
        #extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="none",  # Can be 'none', 'samtools', or 'picard'.
        sort_order="coordinate",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",  # Extra args for samtools/picard sorts.
    threads: 10
    wrapper:
        "v7.6.0/bio/bwa-mem2/mem"


# This rule indexes the combined SCG and TE library using BWA2 for read mapping
rule index_library_for_mapping:
    input:
        "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta"
    output:
        "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta.0123",
        "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta.amb",
        "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta.ann",
        "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta.bwt.2bit.64",
        "{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta.pac",
    log:
        "{species}/processed/dynamics/lib/{feature_library}_and_scg_bwa_index.log"
    message: "Indexing SCG and Feature library {input} with BWA2"
    wrapper:
        "v7.2.0/bio/bwa-mem2/index"


# Rule: Convert SAM to BAM
rule convert_sam_to_bam_reads_to_library:
    # 2 Convert SAM to BAM
    input:
        "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.sam"
    output:
        bam=temp("{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.unsorted.bam") # needs to be called .unsorted.bam otherwise snakemake had problems with unambigous names
    message: "Converting SAM to BAM for {input}"
    threads: 8
    wrapper:
        "v7.5.0/bio/samtools/view"


# Rule: Convert SAM to BAM
rule remove_unmapped:
    # 2 Convert SAM to BAM
    input:
        "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.unsorted.bam"
    output:
        bam=temp("{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.unsorted.no_unmapped.bam") # needs to be called .unsorted.bam otherwise snakemake had problems with unambigous names
    message: "Removing unmapped reads from BAM file for {input}"
    params:
        extra="-b -F 4",  # optional params string
    threads: 2
    wrapper:
        "v7.5.0/bio/samtools/view"

# Rule: Sort BAM file
rule  sort_bam_reads_to_library:
    input:
        "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.unsorted.no_unmapped.bam"
    output:
        "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.sorted.bam"
    message: "Sorting BAM file for {input}"
    log:
        "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg_sort_bam.log",
    threads: 8
    wrapper:
        "v7.5.0/bio/samtools/sort"

# rule deduplicate_bam_with_dedup:
#     input:
#         bam         = "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.sorted.bam"
#     output:
#         dedup_folder= directory("{species}/processed/dynamics/mapped/deduplication/{individual}/"),
#         dedup_bam   = "{species}/processed/dynamics/mapped/deduplication/{individual}/{individual}_scg_feature_library.sorted_rmdup.bam",
#         dedup_hist  = "{species}/processed/dynamics/mapped/deduplication/{individual}/{individual}_scg_feature_library.sorted.hist",
#         dedup_json  = "{species}/processed/dynamics/mapped/deduplication/{individual}/{individual}_scg_feature_library.sorted.dedup.json",
#     message:
#         "Deduplicating BAM file for {input.bam} using dedup for individual {wildcards.individual} in species {wildcards.species}",
#     log: 
#         "{species}/processed/dynamics/mapped/deduplication/{individual}_dedup.log",
#     conda:
#         "../envs/dedup.yaml"
#     resources:
#         mem_mb = 20000   # request 10 GB from cluster / cgroups
#     shell:
#         """
#         mkdir -p {output.dedup_folder}
#         # Set explicit heap size via -Xms (initial) and -Xmx (max)
#         dedup -Xms5g -Xmx20g --input {input.bam} --merged --output {output.dedup_folder} > {log}
#         """

# rule move_deduplicated_to_library_mapped:
#     input:
#         dedup_bam = "{species}/processed/dynamics/mapped/deduplication/{individual}/{individual}_scg_feature_library.sorted_rmdup.bam",
#     output:
#         moved_bam = temp("{species}/processed/dynamics/mapped/{individual}_scg_feature_library.unsorted.dedupped.bam")
#     shell:
#         """
#         mv {input.dedup_bam} {output.moved_bam}
#         """

# # Rule: Sort BAM file
# rule sort_dedup_bam_reads_to_library:
#     # 3 Sort BAM
#     input:
#         "{species}/processed/dynamics/mapped/{individual}_scg_feature_library.unsorted.dedupped.bam"
#     output:
#         "{species}/processed/dynamics/mapped/{individual}_scg_feature_library.sorted.dedupped.bam"
#     message: "Sorting BAM file for {input}"
#     log:
#         "{species}/processed/dynamics/mapped/{individual}.sorted.dedupped.bam.log",
#     threads: 10
#     wrapper:
#         "v7.5.0/bio/samtools/sort"

# Rule: Index BAM file
# SAMTOOLS doesn’t parallelize the indexing work — it only parallelizes compression/decompression.
rule  index_bam_reads_to_library:
    # 4 Index BAM
    input:
        "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.sorted.bam"
    output:
        "{species}/processed/dynamics/{feature_library}/mapped/{individual}_{feature_library}_and_scg.sorted.bam.bai"
    message: "Indexing BAM file for {input}"
    params:
        extra="",  # optional params string
    threads: 5
    wrapper:
        "v7.5.0/bio/samtools/index"

# # Rule: Index BAM file
# # SAMTOOLS doesn’t parallelize the indexing work — it only parallelizes compression/decompression.
# rule  index_dedup_bam_reads_to_library:
#     # 4 Index BAM
#     input:
#         "{species}/processed/dynamics/mapped/{individual}_scg_feature_library.sorted.dedupped.bam" 
#     output:
#         "{species}/processed/dynamics/mapped/{individual}_scg_feature_library.sorted.dedupped.bam.bai"
#     message: "Indexing BAM file for {input}"
#     params: 
#         extra="",  # optional params string
#     threads: 5
#     wrapper:
#         "v7.5.0/bio/samtools/index"