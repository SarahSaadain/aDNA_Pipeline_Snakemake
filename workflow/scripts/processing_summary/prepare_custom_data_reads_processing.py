import json
import pandas as pd

# Snakemake variables are automatically available
reads_file = snakemake.input.reads
dedup_file = snakemake.input.dedup
individual = snakemake.wildcards.individual
output_file = snakemake.output[0]

# Load reads CSV
reads_df = pd.read_csv(reads_file)

# Filter rows for this individual
reads_df = reads_df[reads_df["individual"] == individual]

# Aggregate read-processing steps (sum across samples)
reads_summary = reads_df[[
    "raw_count",
    "adapter_removed_count",
    "quality_filtered_count"
]].sum().to_dict()

# Load DeDup JSON
with open(dedup_file) as f:
    dedup_data = json.load(f)["metrics"]

# Build final table
summary = {
    "individual": individual,
    "raw_reads": reads_summary["raw_count"],
    "after_adapter_removed": reads_summary["adapter_removed_count"],
    "after_quality_filter": reads_summary["quality_filtered_count"],
    "mapped_endogenous_reads": dedup_data["mapped_reads"],
    "endogenous_duplicates_removed": dedup_data["mapped_reads"] - dedup_data["total_removed"]
}

# Write TSV
pd.DataFrame([summary]).to_csv(output_file, sep="\t", index=False)
