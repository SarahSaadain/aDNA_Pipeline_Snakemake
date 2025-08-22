import os
import subprocess
import math
import shutil
import common_aDNA_scripts as common_adna
import ref_genome_processing.common_ref_genome_processing_helpers as common_rgp
import ref_genome_processing.convert_mapped_sam2bam as convert_sam2bam

from multiprocessing import Pool, cpu_count

def execute_mapdamage(sorted_bam_path: str, ref_genome_path: str, processed_output_folder: str, result_output_folder: str, max_threads_per_subprocess: int = 1):

    common_adna.print_debug("Entering execute_mapdamage function")

    pid = os.getpid()
    
    common_adna.print_info(f"[PID {pid}] Running mapDamage2 on {sorted_bam_path} ...")

    if not os.path.exists(sorted_bam_path):
        raise Exception(f"[PID {pid}] BAM file {sorted_bam_path} does not exist!")

    if not os.path.exists(ref_genome_path):
        raise Exception(f"[PID {pid}] Reference genome file {ref_genome_path} does not exist!")

    if os.path.exists(os.path.join(processed_output_folder, "misincorporation.txt")) or os.path.exists(os.path.join(result_output_folder, "misincorporation.txt")):
        common_adna.print_skipping(f"[PID {pid}] mapDamage output already exists in {result_output_folder}!")
    else:

        command = [
            common_adna.PROGRAM_PATH_MAPDAMAGE,
            "-i", sorted_bam_path,
            "-r", ref_genome_path,
            "--folder", processed_output_folder,
            "--merge-reference-sequences",
            "--rescale"
        ]

        try:
            common_adna.run_command(
                command, 
                description=f"[PID {pid}] Running mapDamage on {common_adna.get_filename_from_path(sorted_bam_path)}"
                )
        except subprocess.CalledProcessError:
            return
    
    
    sorted_bam_file_name = common_adna.get_filename_from_path(sorted_bam_path)
    
    # this should be the filename and path created by mapdamage
    rescaled_bam_file_name = sorted_bam_file_name.replace(common_adna.FILE_ENDING_BAM, f".rescaled{common_adna.FILE_ENDING_BAM}")
    rescaled_bam_file_path = os.path.join(processed_output_folder, rescaled_bam_file_name)
    sorted_rescaled_bam_file_path = common_rgp.get_rescaled_bam_path_for_sorted_bam_path(sorted_bam_path)

    if os.path.exists(sorted_rescaled_bam_file_path):
        common_adna.print_skipping(f"[PID {pid}] Rescaled BAM file {rescaled_bam_file_path} already processed!")
        return

    if not os.path.exists(rescaled_bam_file_path):
        common_adna.print_error(f"[PID {pid}] Rescaled BAM file {rescaled_bam_file_path} does not exist!")
        return

    common_adna.print_info(f"[PID {pid}] Sorting and indexing rescaled BAM file {rescaled_bam_file_path} ...")
    convert_sam2bam.execute_sort_bam(rescaled_bam_file_path, sorted_rescaled_bam_file_path, detlete_unsorted_bam=False, threads=max_threads_per_subprocess)

    if not os.path.exists(sorted_rescaled_bam_file_path):
        common_adna.print_error(f"[PID {pid}] Sorted rescaled BAM file {sorted_rescaled_bam_file_path} does not exist!")
        return

    convert_sam2bam.execute_index_bam(sorted_rescaled_bam_file_path, threads=max_threads_per_subprocess)

    sorted_rescaled_bai_path = sorted_rescaled_bam_file_path + common_adna.FILE_ENDING_BAI

    if not os.path.exists(sorted_rescaled_bai_path):
        common_adna.print_error(f"[PID {pid}] Index for rescaled BAM file {sorted_rescaled_bam_file_path} does not exist!")
        return

    common_adna.print_info(f"[PID {pid}] Sorted & rescales BAM file moved to {sorted_rescaled_bam_file_path} ...")
    common_adna.print_info(f"[PID {pid}] Removing unsorted rescaled BAM file {rescaled_bam_file_path} ...")
    os.remove(rescaled_bam_file_path)

    damage_processed_files = common_adna.get_files_in_folder_matching_pattern(processed_output_folder, "*")

    common_adna.print_info(f"[PID {pid}] Moving {len(damage_processed_files)} processed files to {result_output_folder} ...")

    for file in damage_processed_files:
        src_path = file
        dst_path = os.path.join(result_output_folder, common_adna.get_filename_from_path(file))
        if os.path.isfile(src_path):
            shutil.move(src_path, dst_path)

    common_adna.print_info(f"[PID {pid}] Finished running mapDamage on {sorted_bam_path}")
    

def run_mapdamage_for_species(species: str):
    common_adna.print_species_execution(f"Running mapDamage for species {species} ...")

    if not common_adna.config.is_process_step_enabled(
        common_adna.PipelineStages.REFERENCE_GENOME_PROCESSING, 
        common_adna.ReferenceGenomeProcessingSteps.DAMAGE_ANALYSIS,
        species
    ):
        common_adna.print_skipping(f"Creating consensus sequence for species {species} is disabled in the config. Skipping this species.")
        return

    try:
        ref_genome_list = common_rgp.get_reference_genome_file_list_for_species(species)
    except Exception as e:
        common_adna.print_error(f"Failed to get reference genome files for species {species}: {e}")
        return

    for ref_genome_id, ref_genome_path in ref_genome_list:
        common_adna.print_debug(f"Reference genome: {ref_genome_id}")

        mapped_folder = common_adna.get_folder_path_species_processed_refgenome_mapped(species, ref_genome_id)
        list_of_bam_files = common_adna.get_files_in_folder_matching_pattern(mapped_folder, f"*{common_adna.FILE_ENDING_SORTED_BAM}")

        if len(list_of_bam_files) == 0:
            common_adna.print_warning(f"No mapped BAM files found in {mapped_folder} for species {species}.")
            return

        common_adna.print_debug(f"Found {len(list_of_bam_files)} BAM files for species {species}")
        common_adna.print_debug(f"BAM files: {list_of_bam_files}")

        # Calculate the number of threads per subprocess
        # This is used for sorting and indexing the rescaled BAM files
        unused_threads = common_adna.THREADS_DEFAULT - len(list_of_bam_files)
        max_threads_per_subprocess = max(1, math.floor(unused_threads / len(list_of_bam_files)))

        # Prepare parallel task list
        tasks = []
        for sorted_bam_file in list_of_bam_files:
            if not os.path.exists(sorted_bam_file):
                common_adna.print_warning(f"Sorted BAM file {sorted_bam_file} does not exist.")
                continue

            individual_id = common_adna.get_filename_from_path(sorted_bam_file).split(".")[0]
            processed_output_folder = common_adna.get_folder_path_species_processed_refgenome_damage_individual(
                species, ref_genome_id, individual_id)
            result_output_folder = common_adna.get_folder_path_species_results_refgenome_damage_individual(
                species, ref_genome_id, individual_id)
            
            tasks.append((sorted_bam_file, ref_genome_path, processed_output_folder, result_output_folder, max_threads_per_subprocess))

        if tasks:
            num_processes = min(common_adna.THREADS_DEFAULT, cpu_count(), len(tasks))
            common_adna.print_info(f"Running mapDamage with max {num_processes} threads ...")

            with Pool(processes=num_processes) as pool:
                pool.starmap(execute_mapdamage, tasks)

    common_adna.print_success(f"mapDamage analysis for species {species} complete")


def all_species_run_mapdamage():
    common_adna.print_step_execution("Running mapDamage for all species")

    if not common_adna.config.is_process_step_enabled(
        common_adna.PipelineStages.REFERENCE_GENOME_PROCESSING, 
        common_adna.ReferenceGenomeProcessingSteps.DAMAGE_ANALYSIS
    ):
        common_adna.print_skipping("mapDamage analysis is not enabled in the config. Skipping this step.")
        return

    for species in common_adna.FOLDER_SPECIES:
        run_mapdamage_for_species(species)

    common_adna.print_info("mapDamage complete for all species")