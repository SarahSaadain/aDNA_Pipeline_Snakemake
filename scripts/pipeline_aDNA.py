from common_aDNA_scripts import *

#load individual scripts to run within the pipeline
import raw_reads_processing.quality_checking.execute_fastqc as execute_fastqc
import raw_reads_processing.quality_checking.execute_multiqc as execute_multiqc
import raw_reads_processing.execute_fastp_adapter_remove_and_merge as execute_fastp_adapter_remove_and_merge
import raw_reads_processing.polish_fastp_quality_filter as polish_fastp_quality_filter
import raw_reads_processing.polish_fastp_deduplication as polish_fastp_deduplication
import raw_reads_processing.quality_checking.generate_quality_check_report as generate_quality_check_report
import raw_reads_processing.merge_reads_by_individual as merge_reads_by_individual
import raw_reads_processing.analysis.determine_reads_processing_result as determine_reads_processing_result
import raw_reads_processing.analysis.determine_read_length_distribution as determine_read_length_distribution
import raw_reads_processing.analysis.generate_plots_raw_reads_processing as generate_plots_raw_reads_processing
import raw_reads_processing.analysis.contamination.check_contamination_centrifuge as check_contamination_centrifuge
import raw_reads_processing.analysis.contamination.check_contamination_kraken as check_contamination_kraken
import raw_reads_processing.analysis.contamination.check_contamination_ecmsd as check_contamination_ecmsd

import ref_genome_processing.prepare_ref_genome_for_mapping as prepare_ref_genome_for_mapping
import ref_genome_processing.map_aDNA_to_refgenome as map_aDNA_to_refgenome
import ref_genome_processing.create_consensus_sequence_and_mafs as create_consensus_sequence_and_mafs
import ref_genome_processing.analysis.determine_endogenous_reads as determine_endogenous_reads
import ref_genome_processing.analysis.extract_special_sequences as extract_special_sequences
import ref_genome_processing.analysis.determine_coverage_depth_and_breadth as determine_coverage_depth_and_breadth
import ref_genome_processing.analysis.generate_plots_ref_genome_processing as generate_plots_ref_genome_processing
import ref_genome_processing.analyze_and_rescale_damage as analyze_damage

import additional_analysis.species_comparison.analysis.generate_plots_species_compare as generate_plots_species_compare
import additional_analysis.mtdna_analysis.pipeline_mtdna_analysis as pipeline_mtdna_analysis


def run_pipeline_reference_genome_processing():

    print_stage_execution("Starting reference genome processing pipeline ...")

    print_debug(config.get_pipeline_settings(PipelineStages.REFERENCE_GENOME_PROCESSING))

    if not config.is_process_stage_enabled(PipelineStages.REFERENCE_GENOME_PROCESSING):
        print_skipping("Reference genome processing is not enabled in the config. Skipping this stage.")
        return
    
    ############################################################
    # Mapping to reference genome
    ############################################################

    # prepare reference genome for mapping
    prepare_ref_genome_for_mapping.all_species_prepare_ref_genome()
    
    # map reads to reference genome
    map_aDNA_to_refgenome.all_species_map_aDNA_to_refgenome()


    # run mapDamage on the mapped reads
    # this step analyzes the damage patterns in the mapped reads
    analyze_damage.all_species_run_mapdamage()

    # create consensus sequence from mapped reads
    # this step creates a consensus sequence from the mapped reads
    create_consensus_sequence_and_mafs.all_species_create_consensus_sequence_and_mafs()

    # quality control for mapped reads
    # determine endogenous reads
    determine_endogenous_reads.all_species_determine_endogenous_reads()

    # determine coverage depth and breadth
    determine_coverage_depth_and_breadth.all_species_determine_coverage_depth_and_breath()

    # extract special sequences 
    extract_special_sequences.all_species_extract_special_sequences()

    # generate plots for all species to visualize results
    # these contain 
    # 1. coverage depth and breadth
    # 2. endogenous reads
    generate_plots_ref_genome_processing.all_species_generate_plots()


def run_pipeline_raw_reads_processing():

    print_stage_execution("Starting raw reads processing pipeline ...")

    print_debug(config.get_pipeline_settings(PipelineStages.RAW_READS_PROCESSING))

    if not config.is_process_stage_enabled(PipelineStages.RAW_READS_PROCESSING):
        print_skipping("Raw reads processing is not enabled in the config. Skipping this stage.")
        return

    ############################################################
    # Processing of reads
    ############################################################

    # quality control for raw reads using fastqc and multiqc
    execute_fastqc.all_species_fastqc_raw()
    execute_multiqc.all_species_multiqc_raw()

    # adapter removal
    # this step uses fastp to remove adapters from the raw reads
    # it can handle single and paired end reads. paired end reads are merged
    execute_fastp_adapter_remove_and_merge.all_species_fastp_adapter_remove_and_merge()

    # quality control for adapter removed reads using fastqc and multiqc
    execute_fastqc.all_species_fastqc_adapter_removed()
    execute_multiqc.all_species_multiqc_adapter_removed()

    # apply quality filtering to adapter removed reads
    # this step uses fastp to apply quality filtering to the adapter removed reads
    polish_fastp_quality_filter.all_species_fastp_quality_filter()

    # quality control for quality filtered reads using fastqc and multiqc
    execute_fastqc.all_species_fastqc_quality_filtered()
    execute_multiqc.all_species_multiqc_quality_filtered()

    # remove duplicates from quality filtered reads
    polish_fastp_deduplication.all_species_fastp_deduplication()

    # quality control for duplicates removed reads
    execute_fastqc.all_species_fastqc_duplicates_removed()
    execute_multiqc.all_species_multiqc_duplicates_removed()

    # generate quality check report (html) to easily access all qc results
    generate_quality_check_report.all_species_generate_quality_check_report()

    # prepare species for mapping to reference genome
    # this step uses the processed reads and prepares them for mapping to the 
    # reference genome by concatenating the reads and creating different fastq files.
    # fastq files are created per individual
    merge_reads_by_individual.all_species_merge_reads_by_individual()

    # determine reads processing before and after
    determine_reads_processing_result.all_species_determine_determine_reads_processing_result()

    # determine read length distribution
    determine_read_length_distribution.all_species_determine_read_length_distribution()

    # determine contamination using centrifuge
    check_contamination_centrifuge.all_species_run_centrifuge()

    # determine contamination using kraken
    check_contamination_kraken.all_species_run_Kraken()

    # determine contamination using eCMSD
    check_contamination_ecmsd.all_species_run_ecmsd()

     # generate plots for all species to visualize results
    # these contain 
    # 1. reads processing results before and after
    # 2. sequence length distribution
    # 3. contamination analysis
    generate_plots_raw_reads_processing.all_species_generate_plots()

def run_pipeline_post_processing():

    print_stage_execution("Starting post processing pipeline ...")

    print_debug(config.get_pipeline_settings(PipelineStages.POST_PROCESSING))

    if not config.is_process_stage_enabled(PipelineStages.POST_PROCESSING):
        print_skipping("Post processing is not enabled in the config. Skipping this stage.")
        return

    ############################################################
    # Post processing
    ############################################################
    # determine coi
    pipeline_mtdna_analysis.pipeline_mtdna_analysis()

    ############################################################
    # Generate plots
    ############################################################
    # generate comparison plots for species
    # these contain
    # 1. reads processing results before and after
    # 2. coverage depth and breadth
    # 3. endogenous reads
    generate_plots_species_compare.species_generate_comparison_plots()
    

def run_pipeline():

    pid = os.getpid()

    print_stage_execution(f"Starting pipeline (Main PID: {pid}) ...")
    
    ############################################################
    # Processing of reads
    ############################################################
    run_pipeline_raw_reads_processing()    

    ############################################################
    # Mapping to reference genome
    ############################################################
    run_pipeline_reference_genome_processing()

    ############################################################
    # Post processing
    ############################################################
    run_pipeline_post_processing()

    print_success("Pipeline completed successfully.")


def main():
    run_pipeline()

if __name__ == "__main__":  
    main()