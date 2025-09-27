# work in progress


# rule samtools_filter_unmapped:
#     """
#     Apply QC filtration after raw alignment. Three filters are applied concurrently:
#     - (user-defined) minimum base quality
#     - (user-defined) minimum length
#     - (constant)     remove unmapped
#     """
#     input:
#         sam        = assign_aligner_algorithm,
#         reference = rules.index_reference_genome.input.refgen,
#         bwt        = rules.index_reference_genome.output.bwt
#     output:
#         bam        = temp("results/01-preprocess/03-filter/{sample}/{run}/{sample}.bam")
#     params:
#         min_MQ     = config["preprocess"]["filter"]["min-MQ"],
#         min_length = config["preprocess"]["filter"]["min-length"],
#     log:       "logs/01-preprocess/samtools_filter_unmapped/{sample}/{run}.log"
#     benchmark: "benchmarks/01-preprocess/samtools_filter_unmapped/{sample}/{run}.tsv"
#     conda:     "../envs/samtools-1.15.yml"
#     priority:  6
#     threads:   8
#     shell: """
#         samtools view \
#         --threads {threads} \
#         --reference {input.reference} \
#         -q {params.min_MQ} \
#         -e 'length(seq)>{params.min_length}' \
#         -F4 -Sb \
#         {input.sam} > {output.bam} 2> {log}
#     """


rule picard_rmdup:
    """
    Remove PCR duplicates from a BAM file using picard.
    """
    input:
        bams     = "bam",
    output:
        bam     = "results/01-preprocess/06-dedup/picard/{sample}/{sample}.srt.rmdup.bam",
        metrics = "results/01-preprocess/06-dedup/picard/{sample}/{sample}.rmdup.metrics.txt"
    log:       "logs/01-preprocess/picard_rmdup/{sample}.log"
    benchmark: "benchmarks/01-preprocess/picard_rmdup/{sample}.tsv"
    message:    "Removing duplicates with Picard for {input.bams}"
    params:
        extra="--REMOVE_DUPLICATES true --VALIDATION_STRINGENCY LENIENT --ASSUME_SORT_ORDER coordinate",
    wrapper:
        "v7.5.0/bio/picard/markduplicates"


rule samtools_markdup:
    input:
        aln="{sample}.bam",
    output:
        bam="{sample}.markdup.bam",
        idx="{sample}.markdup.bam.csi",
    log:
        "{sample}.markdup.log",
    benchmark: "benchmarks/01-preprocess/samtools_rmdup/{sample}.tsv"
    message: "Removing duplicates with samtools markdup for {input.aln}"
    params:
        extra="-r -s -c --no-PG",
    threads: 2
    wrapper:
        "v7.5.0/bio/samtools/markdup"
