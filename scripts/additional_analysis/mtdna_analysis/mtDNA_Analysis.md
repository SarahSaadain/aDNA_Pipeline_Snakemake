# mtDNA Analysis

The scripts in this folder are used to analyze mtDNA data.

The process is as follows:

1. mtDNA reads are manually provided in the `species/raw/mtdna/` folder as fastq files.
2. The mtDNA reads are aligned to the reference genome using the `scripts/additional_analysis/mtdna_analysis/determine_mtdna_step1_map_to_ref_genome.py` script. The mapped reads are stored in the `species/processed/mtdna/mapped/` folder.
3. The aligned mtDNA reads are used to determine the mtDNA regions using the `scripts/additional_analysis/mtdna_analysis/determine_mtdna_step2_determine_regions.py` script. The region file is stored in the `species/results/mtdna/regions/` folder.
4. A consensus sequence is create for each mapped species using the `scripts/additional_analysis/mtdna_analysis/determine_mtdna_step3_create_and_map_consensus_sequence.py` script. The mapped species reads are taken from the `species/processed/mapped/` folder and the consensus sequence is stored in the `species/results/mtdna/consensus_sequences/` folder.


