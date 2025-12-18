import csv

input_csv = snakemake.input.csv
output_tsv = snakemake.output.tsv
individual = snakemake.params.individual

covered = 0
total = 0

with open(input_csv, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        covered += int(row["covered_bases"])
        total += int(row["total_bases"])

uncovered = total - covered

with open(output_tsv, "w") as out:
    out.write("sample\tcovered_bases\tuncovered_bases\n")
    out.write(f"{individual}\t{covered}\t{uncovered}\n")
