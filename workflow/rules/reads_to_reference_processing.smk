# =================================================================================================
#     Reads to Reference Processing Workflow
# =================================================================================================
# This file coordinates the steps for mapping aDNA reads to a reference and downstream analyses.
# Each include statement brings in rules for a specific processing or analysis step.

# Prepare the reference for mapping (indexing, etc.)
include: "reads_to_reference/processing/prepare_reference_for_mapping.smk"

# Map ancient DNA reads to the reference
include: "reads_to_reference/processing/map_reads_to_reference.smk"

# Deduplicate mapped reads
include: "reads_to_reference/processing/deduplication.smk"

# Analyze DNA damage patterns and rescale BAM files
include: "reads_to_reference/processing/analyze_damage_and_rescale_bam.smk"

# Get the final BAM file
include: "reads_to_reference/processing/get_final_bam.smk"

# Determine coverage depth and breadth statistics
include: "reads_to_reference/analytics/statistics/determine_coverage_depth_breadth.smk"

# Calculate additional mapping statistics using samtools
include: "reads_to_reference/analytics/statistics/determine_samtools_stats.smk"

# Assess the number of endogenous reads
include: "reads_to_reference/analytics/statistics/determine_endogenous_reads.smk"

# Plot endogenous reads statistics
include: "reads_to_reference/plotting/plot_endogenous_reads.smk"

# Plot coverage breadth across the reference 
include: "reads_to_reference/plotting/plot_coverage_breadth.smk"

# Plot coverage depth across the reference
include: "reads_to_reference/plotting/plot_coverage_depth.smk"
# =================================================================================================
# End of reads_to_reference_processing.smk
# =================================================================================================