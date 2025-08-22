import os
from common_aDNA_scripts import *


def get_html_list_of_files(species: str, files: list):
    html_list = ""
    for file in files:

        # make path relative to report folder
        file = file[len(get_folder_path_species_results_qc(species)):]

        # remove leading slash
        file = file[1:]

        html_list += f"<li><a href='{file}'>{file}</a></li>"
    return html_list

def species_generate_quality_check_report(species: str):
    print_info("Generating quality check report for species: %s" % species)

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.GENERATE_QUALITY_CHECK_REPORT,
        species
    ):
        print_skipping(f"Generating quality check report for {species} is not enabled in the config. Skipping this step.")
        return

    # get fastqc html pages
    fastqc_raw_html_files = get_files_in_folder_matching_pattern(get_folder_path_species_results_qc_fastqc_raw(species), f"*{FILE_ENDING_HTML}")
    fastqc_trimmed_html_files = get_files_in_folder_matching_pattern(get_folder_path_species_results_qc_fastqc_adapter_removed(species), f"*{FILE_ENDING_HTML}")
    fastqc_quality_filtered_html_files = get_files_in_folder_matching_pattern(get_folder_path_species_results_qc_fastqc_quality_filtered(species), f"*{FILE_ENDING_HTML}")
    fastqc_duplicates_removed_html_files = get_files_in_folder_matching_pattern(get_folder_path_species_results_qc_fastqc_duplicates_removed(species), f"*{FILE_ENDING_HTML}")

    # get multiqc html pages    
    multiqc_raw_html_files = get_files_in_folder_matching_pattern(get_folder_path_species_results_qc_multiqc_raw(species), f"*{FILE_ENDING_HTML}")
    multiqc_trimmed_html_files = get_files_in_folder_matching_pattern(get_folder_path_species_results_qc_multiqc_adapter_removed(species), f"*{FILE_ENDING_HTML}")
    multiqc_quality_filtered_html_files = get_files_in_folder_matching_pattern(get_folder_path_species_results_qc_multiqc_quality_filtered(species), f"*{FILE_ENDING_HTML}")
    multiqc_duplicates_removed_html_files = get_files_in_folder_matching_pattern(get_folder_path_species_results_qc_multiqc_duplicates_removed(species), f"*{FILE_ENDING_HTML}")

    report_folder = get_folder_path_species_results_qc(species)
    report_file = f"quality_check_report_{species}.html"

    # create report as html file
    with open(os.path.join(report_folder, report_file), "w") as report:
        report.write(f"""
        <html>
            <head>
                <title>Quality Check Report: {species}</title>
            </head>
            <body>
                <h1>Quality Check Report: {species}</h1>
                <h2>MultiQC</h2>
                <h3>Raw Reads</h3>
                <ul>
                    {get_html_list_of_files(species, multiqc_raw_html_files)}
                </ul>
                <h3>Trimmed Reads</h3>
                <ul>
                    {get_html_list_of_files(species, multiqc_trimmed_html_files)}    
                </ul>   
                <h3>Quality Filtered Reads</h3>
                <ul>
                    {get_html_list_of_files(species, multiqc_quality_filtered_html_files)}
                </ul>   
                <h3>Duplicates Removed Reads</h3>
                <ul>
                    {get_html_list_of_files(species, multiqc_duplicates_removed_html_files)}    
                </ul>   
                <h2>FastQC</h2>
                <h3>Raw Reads</h3>
                <ul>
                    {get_html_list_of_files(species, fastqc_raw_html_files)}
                </ul>
                <h3>Trimmed Reads</h3>
                <ul>
                    {get_html_list_of_files(species, fastqc_trimmed_html_files)}
                </ul>   
                <h3>Quality Filtered Reads</h3>
                <ul>
                    {get_html_list_of_files(species, fastqc_quality_filtered_html_files)}
                </ul>   
                <h3>Duplicates Removed Reads</h3>
                <ul>
                    {get_html_list_of_files(species, fastqc_duplicates_removed_html_files)}
                </ul>   
            </body>
        </html>
        """)

    print_info(f"Quality check report generated for species: {species}")



def all_species_generate_quality_check_report():
    print_step_execution("Generating quality check report for all species")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.GENERATE_QUALITY_CHECK_REPORT
    ):
        print_skipping(f"Generating quality check report is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES:
        species_generate_quality_check_report(species)

    print_info("Finished generating quality check report for all species")