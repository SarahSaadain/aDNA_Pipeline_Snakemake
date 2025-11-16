import os

def combine_endogenous_files(individual_files, combined_file_path):
    """Combines individual endogenous read files into a single summary file."""
    if not individual_files:
        raise Exception(f"No individual endogenous read files found to combine.")

    individual_files = sorted(individual_files)

    with open(combined_file_path, "w") as combined_file:
        combined_file.write("Filename,MappedReads,TotalReads,Proportion\n")

        for f in individual_files:
            with open(f, "r") as f_in:
                _ = f_in.readline()  # skip header
                data_line = f_in.readline().strip()
                if data_line:
                    combined_file.write(data_line + "\n")

    print(f"Successfully created combined endogenous reads file: {combined_file_path}")


if __name__ == "__main__":
    combine_endogenous_files(
        snakemake.input,   # list of CSVs
        snakemake.output[0]
    )
