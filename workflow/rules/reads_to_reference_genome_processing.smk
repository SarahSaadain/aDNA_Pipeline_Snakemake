# =================================================================================================
#     Reads to Reference Genome Processing Workflow
# =================================================================================================
# This file coordinates the steps for mapping aDNA reads to a reference genome and downstream analyses.
# Each include statement brings in rules for a specific processing or analysis step.

# Step 1: Prepare the reference genome for mapping (indexing, etc.)
include: "reads_to_reference_genome/processing/prepare_ref_genome_for_mapping.smk"

# Step 2: Map ancient DNA reads to the reference genome
include: "reads_to_reference_genome/processing/map_aDNA_to_refgenome.smk"

# Step 3: Analyze DNA damage patterns and rescale BAM files
include: "reads_to_reference_genome/processing/analyze_damage_and_rescale_bam.smk"

# Step 4: Create consensus sequences and calculate minor allele frequencies (MAFs)
include: "reads_to_reference_genome/processing/create_consensus_sequence_and_mafs.smk"

# Step 5: Determine coverage depth and breadth statistics
include: "reads_to_reference_genome/analytics/statistics/determine_coverage_depth_breadth.smk"

# Step 6: Calculate additional mapping statistics using samtools
include: "reads_to_reference_genome/analytics/statistics/determine_samtools_stats.smk"

# Step 7: Assess the number of endogenous reads
include: "reads_to_reference_genome/analytics/statistics/determine_endogenous_reads.smk"

# Step 8: Plot endogenous reads statistics
include: "reads_to_reference_genome/plotting/plot_endogenous_reads.smk"

# Step 9: Plot coverage breadth across the reference genome
include: "reads_to_reference_genome/plotting/plot_coverage_breadth.smk"

# Step 10: Plot coverage depth across the reference genome
include: "reads_to_reference_genome/plotting/plot_coverage_depth.smk"
# =================================================================================================
# End of reads_to_reference_genome_processing.smk
# =================================================================================================