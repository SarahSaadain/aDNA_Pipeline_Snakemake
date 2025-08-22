import os
import subprocess

from common_aDNA_scripts import *

def execute_convert_sam_to_bam(sam_file: str, bam_file: str, sorted_bam: str, threads: int=THREADS_DEFAULT, delete_sam: bool=True, detlete_unsorted_bam: bool=True):

    if not os.path.exists(sorted_bam):

        if not os.path.exists(bam_file):
            print_info(f"Converting {sam_file} to BAM...")

            if not os.path.exists(sam_file):
                print_error(f"Input SAM file {sam_file} does not exist!")
                return

            try:
                command_sam_to_bam = [
                    PROGRAM_PATH_SAMTOOLS, 
                    PROGRAM_PATH_SAMTOOLS_VIEW, 
                    "-@", str(threads), 
                    "-bS", sam_file,
                    "-o", bam_file ]
                
                run_command(
                    command_sam_to_bam, 
                    description=f"Converting {sam_file} to {bam_file}"
                )

                if delete_sam and os.path.exists(sam_file):
                    print_info(f"Removing SAM file {sam_file}...")
                    try:
                        os.remove(sam_file)
                    except Exception as e:
                        print_warning(f"Failed to remove SAM file {sam_file}: {e}")
            except Exception as e:
                print_error(f"Failed to convert {sam_file} to BAM: {e}")
                return
        else:
            print_info(f"Conversion for {sam_file} already exists. Skipping.")

    execute_sort_bam(bam_file, sorted_bam, detlete_unsorted_bam, threads)
    
    execute_index_bam(sorted_bam, threads)
    
    
    print_success(f"Conversion and indexing of {sam_file} completed successfully.")

def execute_sort_bam(source_unsorted_bam_path: str, target_sorted_bam_path: str, detlete_unsorted_bam: bool=True, threads: int=THREADS_DEFAULT):

    print_info(f"Sorting {source_unsorted_bam_path}...")

    if not os.path.exists(source_unsorted_bam_path):
        print_error(f"Input BAM file {source_unsorted_bam_path} does not exist!")
        return
    
    if os.path.exists(target_sorted_bam_path):
        print_skipping(f"Sorting for {source_unsorted_bam_path} already exists. Skipping.")
        return
        
    command_sort = [
        PROGRAM_PATH_SAMTOOLS,
        PROGRAM_PATH_SAMTOOLS_SORT,
        "-@", str(threads), 
        source_unsorted_bam_path, 
        "-o", target_sorted_bam_path ]
    
    try:

        run_command(
            command_sort, 
            description=f"Sorting {source_unsorted_bam_path}",
        )

        print_info(f"Sorting of {source_unsorted_bam_path} completed successfully.")
        
        # Optional cleanup of intermediate BAM file if the sorted was created
        if detlete_unsorted_bam and os.path.exists(target_sorted_bam_path):
            print_info(f"Removing unsorted BAM file {source_unsorted_bam_path}...")
            os.remove(source_unsorted_bam_path)

    except Exception as e:
        print_error(f"Failed to sort {source_unsorted_bam_path}: {e}")
        return
    
def execute_index_bam(bam_path: str, threads: int=THREADS_DEFAULT):

    print_info(f"Indexing {bam_path}...")

    if not os.path.exists(bam_path):
        print_error(f"Input BAM file {bam_path} does not exist!")
        return
    
    # Index the sorted BAM file
    indexed_bam = bam_path + FILE_ENDING_BAI

    if os.path.exists(indexed_bam):
        print_skipping(f"Index for {bam_path} already exists. Skipping.")
    
    try:
        command_index = [
            PROGRAM_PATH_SAMTOOLS, 
            PROGRAM_PATH_SAMTOOLS_INDEX,
            "-@", str(threads),
            bam_path,
            indexed_bam ]

        run_command(
            command_index, 
            description=f"Indexing {bam_path}"
        )
        
        print_info(f"Indexing of {bam_path} completed successfully.")
    except Exception:
        return
        