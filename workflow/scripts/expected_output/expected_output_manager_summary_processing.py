
def get_expected_output_summary_processing(species):
    expected_output = []

    expected_output.append(os.path.join(species, "results", "summary", f"{species}_multiqc.overall.html"))

    #{species}/results/analytics/{individual}_multiqc.html
    for individual in get_individuals_for_species(species):
        expected_output.append(os.path.join(species, "results", "summary", f"{individual}_multiqc.html"))

    return expected_output