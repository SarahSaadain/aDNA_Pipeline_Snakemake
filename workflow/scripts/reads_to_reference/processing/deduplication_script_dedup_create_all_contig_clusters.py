import os
import math
import logging

input = snakemake.input
output = snakemake.output
params = snakemake.params

# Read all contigs
with open(input.contigs) as f:
    contigs = [line.strip() for line in f]

n = params.contigs_per_cluster
os.makedirs(output.cluster_folder, exist_ok=True)

# Compute number of groups
num_groups = math.ceil(len(contigs) / n)

# Write one bed file per cluster
# use 1000000000 as the upper bound to ensure full contig coverage. 
# Samtools view will clip to the actual contig length.
for i in range(num_groups):
    start = i * n
    end = min((i + 1) * n, len(contigs))
    group_contigs = contigs[start:end]
    group_file = os.path.join(output.cluster_folder, f"cluster_{start+1}_{end}.bed")  # 1-based group numbering
    with open(group_file, "w") as out:
        logger.info(f"Writing contig cluster file: {group_file} with contigs {start+1} to {end}")
        out.write("\n".join(group_contigs) + "\n")