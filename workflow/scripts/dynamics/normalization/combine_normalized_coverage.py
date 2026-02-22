#!/usr/bin/env python3
import os
import glob
import csv
from snakemake.script import snakemake

def get_individual_name(filename):
    # Extract the individual name from the filename (before _normalized.tsv)
    base = os.path.basename(filename)
    return base.split('_coverage.tsv')[0]

# Find all _normalized.tsv files (recursively)
normalized_files = snakemake.input.coverage_files
output_file = snakemake.output.combined

fle_coverage = {}  # {te_name: {individual: normalized_coverage}}
individuals = []

for file in normalized_files:
    individual = get_individual_name(file)
    individuals.append(individual)
    with open(file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            library_read_type = row['type']
            library_read = row['name']
            coverage = row['normalized_coverage']
            if library_read not in fle_coverage:
                fle_coverage[library_read] = {}
            fle_coverage[library_read][individual] = coverage
            fle_coverage[library_read]["type"] = library_read_type

# Sort individuals and TEs
individuals = sorted(set(individuals))
fle_names = sorted(fle_coverage.keys())

# Write combined output

with open(output_file, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['library_read'] + individuals)
    for library_read in fle_names:
        row = [ f"{library_read}"]
        for ind in individuals:
            row.append(fle_coverage[library_read].get(ind, 'NA'))
        writer.writerow(row)

print(f"Combined normalized coverage written to {output_file}")
