# =================================================================================================
#     Raw Reads Processing Workflow
# =================================================================================================
# This file orchestrates the main steps for processing raw sequencing reads in the aDNA pipeline.
# Each include statement brings in a specific rule or set of rules for a processing or analysis step.

# Step 1: Remove sequencing adapters from raw reads
include: "raw_read/processing/adapter_removal.smk"

# Step 2: Filter reads based on quality thresholds
include: "raw_read/processing/quality_filtering.smk"

# Step 3: Merge reads by individual sample for downstream analysis
include: "raw_read/processing/merge_by_individual.smk"

# Step 4: Generate statistics on processed reads
include: "raw_read/analytics/statistics/reads_processing.smk"

# Step 5: Perform FastQC quality checks on raw and processed reads
include: "raw_read/analytics/quality/fastqc_check.smk"

# Step 6: Aggregate FastQC results using MultiQC
include: "raw_read/analytics/quality/multiqc_check.smk"

# Step 7: Generate a comprehensive quality check report
include: "raw_read/analytics/quality/generate_quality_check_report.smk"

# Step 8: Plot comparison of reads before and after processing
include: "raw_read/plotting/plot_comparison_reads_before_after_processing.smk"

# Step 9: Check for contamination using ECMSD
include: "raw_read/analytics/contamination/check_contamination_ecmsd.smk"

# Step 10: Check for contamination using Centrifuge
include: "raw_read/analytics/contamination/check_contamination_centrifuge.smk"

# Step 11: Check for contamination using Kraken
include: "raw_read/analytics/contamination/check_contamination_kraken.smk"
# =================================================================================================
# End of raw_reads_processing.smk
# =================================================================================================
