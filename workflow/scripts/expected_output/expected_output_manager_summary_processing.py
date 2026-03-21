def get_expected_output_summary_processing(species):
    expected_output = []

    expected_output.append(f"{species}/results/summary/species_level/{species}_multiqc.overall.html")

    for individual in get_individuals_for_species(species):
        expected_output.append(f"{species}/results/summary/individual_level/{individual}_multiqc.html")

    return expected_output