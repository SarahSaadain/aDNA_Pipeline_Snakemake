import os

def determine_endogenous_reads_from_stats(stats_file, output_file):
    """Parse samtools stats file and compute mapped/unmapped proportion."""
    total_reads = 0
    mapped_reads = 0

    with open(stats_file) as f:
        for line in f:
            if line.startswith("SN"):
                parts = line.strip().split("\t")
                if parts[1].startswith("raw total sequences:"):
                    total_reads = int(parts[2])
                elif parts[1].startswith("reads mapped:"):
                    mapped_reads = int(parts[2])

    proportion = mapped_reads / total_reads if total_reads > 0 else 0.0

    with open(output_file, "w") as out:
        out.write("Filename,MappedReads,TotalReads,Proportion\n")
        out.write(f"{os.path.basename(stats_file)},{mapped_reads},{total_reads},{proportion:.4f}\n")


if __name__ == "__main__":
    determine_endogenous_reads_from_stats(
        snakemake.input.stats,   # input stats file
        snakemake.output.csv      # output CSV
    )
