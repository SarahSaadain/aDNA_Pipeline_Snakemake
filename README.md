# aDNA Pipeline

This project contains a pipeline to analyze raw ancient data, obtained from the sequencing facility. The pipeline includes various Snakemake workflows to process, analyze, and generate reports on the sequence quality, which helps decide if an aDNA extraction and sequencing was successful, and further polishes the data for downstream analyses.

Note: This pipeline is still in development, but can already be used for analysis.

## Setup Overview

- Before running the pipeline, ensure you have an environment with Snakemake and the required dependencies installed.
- Required dependencies for pipeline processing will be installed automatically, except for the contamination analysis tools. Those have to be installed separately and their details need to be added to the config file.
  - The pipeline supports ECMSD for contamination analysis. Ensure ECMSD is configured in the `config.yaml` file under `contamination_analysis`.
  - The pipeline supports Centrifuge for contamination analysis. Ensure Centrifuge is configured in the `config.yaml` file under `contamination_analysis`.
- You need to add species details to the pipeline.

When adding a new species, make sure to 
- the species folder should be placed in the root folder of your pipeline
- add the folder name should match the species key which is defined in `config.yaml` below `species:` 
- your reads should be renamed according to the naming convention specified in the [RAW Reads filenames](#RAW-Reads-filenames) section
- the pipeline supports automatically moving the raw reads to the `<species>/raw/reads/` folder as well as the reference to the `<species>/raw/ref/` folder. Simply provide the files in the `<species>` folder. Alternatively, you can manually move the files to the respective folders.
  - provide the raw reads in `<species>/raw/reads/` folder
  - provide the reference in `<species>/raw/ref/` folder
- all other folders will be created and populated automatically
  - folder `<species>/processed/` contains the intermediary files during processing. Most of these files are marked as temporary and will be deleted at the end of the pipeline. Some files are kept to allow reprocessing the pipeline from different points in case something fails.
  - folder `<species>/results/` contains the final results and reports. 
  - general reads processing data will be in either `processed`or `results`. Everything related to a reference will have a `<reference>` folder under `processed`or `results`. Typically, only the `results` folder will contain information required for further analyis. In case more information is required, the original files can often be found in the `processed` folder. Some exemptions include `*.sam` and unsorted `*.bam` files. These are deleted to save storrage space. Most other files are kept in order to allow reprocessing the pipeline from different points in case something fails. If a step should be repeated, the relevant files need to be deleted manually. 

## Configuration File Structure for aDNA Pipeline (`config.yaml`)

The `config.yaml` file is used to configure the aDNA pipeline. It contains settings such as project name, the species list and the pipeline stages and process steps.

### Global Settings

* **project\_name**: Name of the aDNA project.

### Pipeline Settings

Defines the overall pipeline behavior, including execution controls and process details.

#### Pipeline Stages and Process Steps

* The pipeline is broken into **stages** (e.g., `raw_reads_processing`, `reference_processing`, `post_processing`).
* Each stage contains multiple **process steps** (e.g., `adapter_removal`, `deduplication`, ...).
* Both stages and process steps can be controlled with `execute: true/false` flags to enable or disable them.
* Some process steps include additional configurable settings (e.g., adapter sequences, database paths, ...).

#### Important Defaults

* You **do not need to specify all stages or process steps** explicitly.
* Any **stage or process step not provided in the config defaults to `execute: true`** and will be executed.

### Example `config.yaml`

```yaml
# config.yaml - Configuration file for aDNA pipeline
# This file contains settings for various stages of the pipeline


project_name: "aDNA_Project"

# Pipeline stages and their configurations
pipeline:

  # Global settings
  global:
    # When true, existing output files will be skipped to avoid re-computation (Default: true)
    skip_existing_files: true

  # Stages of the pipeline

  # Raw reads processing
  # Includes quality checking, adapter removal, quality filtering, merging, 
  # contamination analysis, and statistical analysis
  raw_reads_processing:
    # When true, this stage will be executed. (Default: true)
    execute: true

    # Sub-stages with their respective settings
    # Quality checking of raw reads
    quality_checking_raw:
      # When true, this sub-stage will be executed (Default: true)
      execute: true
    
    # Adapter removal from raw reads
    adapter_removal:
      # When true, this sub-stage will be executed (Default: true)
      execute: true

      # Settings for adapter removal
      settings: 
        # Minimum quality score for adapter removal
        min_quality: 0
        # Minimum length of reads after adapter removal
        min_length: 0
        # Optional: Adapter sequences for read 1 and read 2
        # If not provided, fastp will try to identify adapters automatically
        adapters_sequences:
          # Adapter sequence for read 1
          r1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
          # Adapter sequence for read 2
          r2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" 
    
    # Quality checking of trimmed reads
    quality_checking_trimmed:
      # When true, this sub-stage will be executed (Default: true)
      execute: true

    # Quality filtering of trimmed reads
    quality_filtering:
      # When true, this sub-stage will be executed (Default: true)
      execute: true
      settings:
        # Minimum quality score for quality filtering
        min_quality: 15
        # Minimum length of reads after quality filtering
        min_length: 30

    # Quality checking of quality-filtered reads
    quality_checking_quality_filtered:
      # When true, this sub-stage will be executed (Default: true)
      execute: true

    # Quality checking of merged reads
    quality_checking_merged:
      # When true, this sub-stage will be executed (Default: true)
      execute: true

    # Contamination analysis
    contamination_analysis:
      # When true, this sub-stage will be executed (Default: true)
      execute: true
      tools:
        # ECMSD tool settings for contamination analysis
        ecmsd:
          # When true, this tool will be executed (Default: true)
          execute: true
          settings:
            # Optional: Path to the conda environment for ECMSD
            # If not provided, the default environment will be used
            #conda_env: "../../../../envs/ecmsd.yaml"
            # Path to the ECMSD executable
            # Curretnly, ecmsd can not be installed via conda. Provide the path to the shell script to run ECMSD.
            executable: "/path/to/ecmsd/shell/ECMSD.sh"
        # Centrifuge tool settings for contamination analysis
        centrifuge:
          # When true, this tool will be executed (Default: true)
          execute: true
          settings:
            # Optional: Path to the conda environment for Centrifuge
            # If not provided, the default environment will be used
            #conda_env: "../../../../envs/centrifuge.yaml"
            # Path to the Centrifuge index
            index: "/path/to/centrifuge_index"
    
    # Statistical analysis
    statistical_analysis:
      # When true, this sub-stage will be executed (Default: true)
      execute: true

  # Reference processing
  reference_processing:
    # When true, this stage will be executed (Default: true)
    execute: true
    
    # Deduplication settings
    deduplication:
      # When true, this sub-stage will be executed (Default: true)
      execute: false
      settings:
        # To increase performance, deduplication will be done per cluster of contigs
        # Below settings define how the contigs will be clustered
        # Optional: Maximum number of contigs per cluster (Default: 10 if not specified)
        max_contigs_per_cluster: 10
        # Optional: Maximum number of contigs per cluster (Default: 500 if not specified)
        max_contigs_per_cluster: 500
    
    # Damage rescaling settings for mapDamage2
    damage_rescaling:
      # When true, this sub-stage will be executed (Default: true)
      execute: true

    # Damage analysis settings for mapDamage2
    damage_analysis:
      # When true, this sub-stage will be executed (Default: true)
      execute: true

    # Endogenous reads analysis settings
    endogenous_reads_analysis: 
      # When true, this sub-stage will be executed (Default: true)
      execute: true

    # Coverage analysis settings
    coverage_analysis: 
      # When true, this sub-stage will be executed (Default: true)
      execute: true

# Species details
species:
  Bger:
    name: "Blatella germanica"
  Dsim:
    name: "Drosophila simulans"
```

## Running the Pipeline

The aDNA pipeline is implemented using Snakemake, a workflow management system. Snakemake ensures reproducibility and efficient execution of the pipeline.

### Running the Pipeline

To run the pipeline, navigate to the root directory containing the `workflow/Snakefile` and execute:

```bash
#snakemake --cores <number_of_threads> --use-conda
snakemake --cores <number_of_threads> --use-conda --keep-going
```

Replace `<number_of_threads>` with the number of CPU threads you want to allocate for the pipeline.

**Note:** 
* The `--use-conda` flag enables the use of conda environments specified in the `Snakefile`.
* The `--keep-going` flag allows the pipeline to continue even if a rule fails. Somtimes the analysis of ECMSD fails due to issues with the input data (e.g. low quality reads or low coverage). In this case, the rest of the pipeline can still be executed.
* The number of threads can be adjusted using the `--cores` option when running Snakemake.

### Running the Pipeline in the Background

Depending on the size of the data, it may take some time to complete the pipeline. Thus it is recommended to run the pipeline in the background. You can do this by running the following command:

```bash
nohup snakemake --cores <number_of_threads> --use-conda --keep-going > pipeline.log 2>&1 &
```

### Restarting the Pipeline

Snakemake automatically tracks the state of the pipeline and will only re-run steps that are incomplete or outdated. If you want to restart the pipeline from the beginning, you can delete the relevant output files and re-run the pipeline.

## Folder Structure

### Species Folders

The project contains folders for different species, which contain the raw data, processed data, and results for each species.

#### RAW Reads Filenames

The pipeline expects input read files to follow a standardized naming convention:

```bash
<Individual>_<Original_Filename>.fastq.gz
```

Following this convention ensures proper organization and automated processing within the pipeline.  

##### Filename Components:
- **`<Individual>`** – A unique identifier for the sample or individual.  
- **`<Original_Filename>`** – The original filename assigned by the sequencing platform.  
- **`.fastq.gz`** – The expected file extension, indicating compressed FASTQ format.  

Notes:
 - This name must contain `_R1_` and, if paired-end, `_R2_`.
 - For paired-end data, the name of the reads must be identical except for `_R1_` and `_R2_`.
 - Individual names/IDs will be used to name the output files as well as in the reports and plots.

#### Example:
```
Bger1_326862_S37_R1_001.fastq.gz
```