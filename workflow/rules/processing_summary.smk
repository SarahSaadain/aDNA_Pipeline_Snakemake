include: "processing_summary/prepare_custom_data_breadth.smk"
include: "processing_summary/prepare_custom_data_depth.smk"
include: "processing_summary/prepare_custom_data_ecmsd.smk"
include: "processing_summary/prepare_custom_data_reads_processing.smk"
include: "processing_summary/prepare_custom_folder_links.smk"


# create_multiqc_species.smk
# Contains rules for generating MultiQC reports for each species.
include: "processing_summary/create_multiqc_species.smk"


# create_multiqc_species_individual.smk
# Contains rules for generating MultiQC reports for each species and individual.
include: "processing_summary/create_multiqc_species_individual.smk"

