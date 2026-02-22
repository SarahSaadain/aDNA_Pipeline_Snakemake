import os
import logging

# -----------------------------------------------------------------------------------------------
# Get all expected output file paths for reference processing
def get_expected_output_reference_processing(species):

    if config.get("pipeline", {}).get("reference_processing", {}).get("execute", True) == False:
        logging.info(f"Skipping reference processing for {species}. Disabled in config.")
        return []

    expected_outputs = []

    try:
    # Get all reference for the species
        references_list = get_references_ids_for_species(species)

    except Exception as e: 
        # Print error if reference files are missing or inaccessible
        logging.error(e)
        return []
    
     # Get all individuals for the species
    individuals = get_individuals_for_species(species)

    for reference in references_list:

        if config.get("pipeline", {}).get("reference_processing", {}).get("endogenous_reads_analysis", {}).get("execute", True) == True:
            # Add endogenous and coverage plots for each reference and individual
            expected_outputs.append(os.path.join(species, "results" ,reference, "plots", "endogenous_reads", f"{species}_{reference}_endogenous_reads_bar_chart.png"))
            expected_outputs.append(os.path.join(species, "results" ,reference, "plots", "endogenous_reads", f"{species}_{reference}_raw_and_endogenous_reads_bar_chart.png"))
        else:
            logging.info(f"Skipping endogenous reads analysis for species {species} and reference {reference}. Disabled in config.")
        
        if config.get("pipeline", {}).get("reference_processing", {}).get("coverage_analysis", {}).get("execute", True) == True:
            expected_outputs.append(os.path.join(species, "results" ,reference, "plots", "coverage", f"{species}_{reference}_individual_depth_coverage_violin.png"))
            expected_outputs.append(os.path.join(species, "results" ,reference, "plots", "coverage", f"{species}_{reference}_individual_depth_coverage_bar.png"))
            expected_outputs.append(os.path.join(species, "results" ,reference, "plots", "coverage", f"{species}_{reference}_individual_coverage_breadth_bar.png"))
            expected_outputs.append(os.path.join(species, "results" ,reference, "plots", "coverage", f"{species}_{reference}_individual_coverage_breadth_violin.png"))

        else:
            logging.info(f"Skipping coverage analysis for species {species} and reference {reference}. Disabled in config.")

        #expected_outputs.append(os.path.join(species, "results" ,reference_id, "multiqc.html"))


        for individual in individuals:
            
            expected_outputs.append(os.path.join(species, "processed" ,reference, "mapped", f"{individual}_{reference}_final.bam"))
            expected_outputs.append(os.path.join(species, "processed" ,reference, "mapped", f"{individual}_{reference}_final.bam.bai"))

            #"{species}/results/{reference}/analytics/{individual}_multiqc.html",
            expected_outputs.append(os.path.join(species, "results" ,reference, "analytics", f"{individual}_{reference}_multiqc.html"))

            if config.get("pipeline", {}).get("reference_processing", {}).get("damage_analysis", {}).get("execute", True) == True:
                #{species}/results/{reference}/analytics/{individual}/mapdamage/misincorporation.txt
                expected_outputs.append(os.path.join(species, "results" ,reference, "analytics", individual, "mapdamage", "misincorporation.txt"))
                #expected_outputs.append(directory(os.path.join(species, "results" ,reference_id, "damage", "damageprofile", ind)))
            else:
                logging.info(f"Skipping damage analysis for species {species} and individual {individual} to reference {reference}. Disabled in config.")

    return expected_outputs
