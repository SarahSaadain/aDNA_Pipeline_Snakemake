# aDNA Pipeline

This project contains a pipeline to analyze raw ancient data, obtained from the sequencing facility. The pipeline includes various scripts to process, analyze, and generate reports on the sequence quality, which helps decide if an aDNA extraction and sequencing was successfull, and further polishes the data for downstream analyses.

Note: This pipeline is still in development, but can already be used for analysis

## Setup Overview

- Before running the pipeline, ensure you have the necessary dependencies installed. Please refer to the [Requirements](#Requirements) section for the necessary dependencies and installation instructions.
- Raw reads and reference genome must be provided in the relevant folders. 
    - Raw reads shoud be renamed according to the naming convention specified in the [RAW Reads filenames](#RAW-Reads-filenames) section. Also see the [Manually renaming the raw reads files](#Manually-renaming-the-raw-reads-files) section.
    - Reference genome must be provided in the `species/raw/ref_genome/` folder.
    - Species must be added to the config file
    - Please refer to the [Species Folders](#Species-Folders) section for the expected folder structure.
- You need to add the project path to the config file.

## Configuration File Structure for aDNA Pipeline (`config.yaml`)

The `config.yaml` file is used to configure the aDNA pipeline. It contains settings such as project name, project description, default number of threads, adapter sequences for adapter removal, species-specific settings, and paths to external tools.

### Global Settings

* **project\_name**: Name of the aDNA project.
* **project\_description**: Brief description of the project.
* **path\_adna\_project**: Path to the main project directory.

### Pipeline Settings

Defines the overall pipeline behavior, including execution controls and process details.

* **threads**: Number of CPU threads to use.
* **log\_level**: Logging verbosity (e.g., INFO, DEBUG).

#### Pipeline Stages and Process Steps

* The pipeline is broken into **stages** (e.g., `raw_reads_processing`, `reference_genome_processing`, `post_processing`).
* Each stage contains multiple **process steps** (e.g., `adapter_removal`, `deduplication`).
* Both stages and process steps can be controlled with `execute: true/false` flags to enable or disable them.
* Some process steps include additional configurable settings (e.g., adapter sequences, database paths).

#### Important Defaults

* You **do not need to specify all stages or process steps** explicitly.
* Any **stage or process step not provided in the config defaults to `execute: true`** and will be executed.
* This applies both to the **global pipeline entry** and to any **species-specific pipeline overrides** (see below).

### Species-specific Settings

Defines species included in the project, each with:

* A unique code (e.g., `Bger`)
* Full species name
* Folder name used for species-specific data

### Species-specific Pipeline Overrides

* Each species can optionally include a **pipeline** section to override the global pipeline settings.
* These overrides allow you to **customize process steps within stages** for that species.
* You **cannot enable or disable entire stages** here—only individual process steps.

### Cross-species Comparisons

* Configures comparative analyses between multiple species.
* Each comparison has a unique name and lists species with their reference genome paths.

Structure:


### Comparison of species results
*   `compare_species`
    * `comparison ID/Name`:  First comparison. Unique Name used for file names
        * `species ID 1`: first species 
            * `species_id`: Optional species id. If not provided, the config species ID will be used instead
            * `reference_genome`: Name of reference genome fo comparison. e.g.: "refgenome.fna"
        * `species ID 2`: second species 
            * `species_id`: Optional species id. If not provided, the config species ID will be used instead
            * `reference_genome`: Name of reference genome fo comparison. e.g.: "refgenome.fna"
        * `species ID ...`: ...
            * `reference_genome`: ...
        * `species ID N`: n-th species 
            * `reference_genome`: Name of reference genome fo comparison. e.g.: "refgenome.fna"
    * `comparison ID/Name`:  Second comparison. Unique Name used for file names
        * `species ID 1`: first species 
        * ...

Note: the filed `species_id` can be used optionally to specify a specific species. This is usefull if you want to compare the same species with different reference genomes. If you want to compare different species, then this value can be skipped if the species ID is provided as the parent config node.


### External Tools

* Defines command names or paths to software tools required by the pipeline (e.g., `fastp`, `bwa`, `samtools`, `kraken`).
* Enables flexibility for different environments or tool versions.

**Note**: If the tools are not provided, default values are used and it is expected that the tool can be called via the command line directly.

### Example `config.yaml`

```yaml
# config.yaml - Configuration file for aDNA pipeline

# Global settings
project_name: "aDNA_Project"
project_description: "Analysis of ancient DNA data"
path_adna_project: "/path/to/project" #Main project path

pipeline:
  threads: 50  # Number of threads to use for the pipeline
  log_level: "INFO" # log level. e.g.: INFO, DEBUG, ERROR, ...

  raw_reads_processing:     # stage
    execute: true           # disable or enable this stage
    quality_checking_raw:   # process step
      execute: true         # disable or enable this process step
    adapter_removal:
      execute: true
      settings:             
        adapters_sequences:
          r1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
          r2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" 
    quality_checking_adapter_removed:
      execute: true
    quality_filtering:
      execute: true
    quality_checking_quality_filtered:
      execute: true
    deduplication:
      execute: true
    quality_checking_deduplication:
      execute: true
    generate_quality_check_report:
      execute: true
    merge_reads_by_individual:
      execute: true
    reads_processing_analysis:
      execute: true
    read_length_distribution_analysis: 
      execute: true
    contamination_analysis:
      execute: true
      tools:
        centrifuge:
          execute: true
          settings:
            database: "/path/to/centrifuge_db"
        kraken: 
          execute: true
          settings:
            database: "/path/to/kraken_db"
        ecmsdb:
          execute: true
    generate_raw_reads_plots: 
      execute: true

  reference_genome_processing:
    execute: true
    prepare_reference_genome:
      execute: true
    map_reads_to_reference_genome:
      execute: true
    create_consensus_sequence:
      execute: true
    damage_analysis:
      execute: true
    endogenous_reads_analysis: 
      execute: true
    coverage_analysis: 
      execute: true
    generate_reference_genome_plots: 
      execute: true

  post_processing:
    execute: true
    mtdna_analysis:
      execute: true
    comparison_plots:
      execute: true

# Species-specific settings
species:
  Bger:
    name: "Blatella germanica"
    folder_name: "Bger"
  Dsim:
    name: "Drosophila simulans"
    folder_name: "Dsim"
  Phortica:
    name: "Phortica"
    folder_name: "Phortica"

# cross-species comparison of depth/breadth, endogenous reads, ...
compare_species:
  Bger_Dsim_comparison: # Unique Name, choose whatever fits. This will be used in the filenames
    Bger: # must be ID from above
      reference_genome: "refgenome.fna"
    Dsim:
      reference_genome: "refgenome.fna"
  Bger_Dsim_Phortica_comparison: # Unique Name
    Bger:
      reference_genome: "refgenome.fna"
    Dsim:
      reference_genome: "refgenome.fna"
    Phortica:
      reference_genome: "refgenome.fna"

# Paths to external tools
tools:
  fastp: "fastp"
  sga: "sga"
  multiqc: "multiqc"
  fastqc: "fastqc"
  bwa: "bwa"
  bedtools: "bedtools"
  samtools: "samtools"
  angsd: "angsd"
  seqkit: "seqkit"
  kraken: "kraken"
  mapdamage: "mapDamage"
```

## Folder Structure

### Species Folders

The project contains folders for different species, which contain the raw data, processed data, and results for each species.

When adding a new species, make sure to 
- add the folder name to the `config.yaml`
- provide the raw reads in `<species>/raw/reads/` folder
- provide the reference genome in `<species>/raw/ref_genome/` folder
- provide mtDNA reads in `<species>/raw/mtdna/` folder
- all other folders will be created and populated automatically
  - folder `<species>/processed/` contains the intermediary files during processing
  - folder `<species>/results/` contains the final results and reports
  - general reads processing data will be in either `processed`or `results`. Everything related to a reference genome will have a `<reference_genome>` folder under `processed`or `results`. Typically, only the `results` folder will contain information required for further analyis. In case more information is required, the original files can often be found in the `processed` folder. Some exemptions include `*.sam` and unsorted `*.bam` files. These are deleted to save storrage space. Most other files are kept in order to allow reprocessing the pipeline from different points in case something fails. If a step should be repeated, the relevant files need to be deleted manually. 

#### RAW Reads Filenames

The pipeline expects input read files to follow a standardized naming convention:

```bash
<Individual>_<Protocol>_<Original_Filename>.fastq.gz
```

Following this convention ensures proper organization and automated processing within the pipeline.  

##### Filename Components:
- **`<Individual>`** – A unique identifier for the sample or individual.  
- **`<Protocol>`** – The sequencing or library prep protocol (to enable easy comparison of different protocols or sequencing technologies) (optional).  
- **`<Original_Filename>`** – The original filename assigned by the sequencing platform.  
- **`.fastq.gz`** – The expected file extension, indicating compressed FASTQ format.  

#### Example:

```
Bger1_S_326862_S37_R1_001.fastq.gz
```

#### Manually renaming the raw reads files

For manually renaming the raw read files, use the `rename.py` script located at https://github.com/SarahSaadain/aDNA_Project_Data/blob/main/resources/rename.py .  

### Scripts Folder

The `scripts/` folder contains all necessary scripts for the aDNA pipeline, organized into subfolders corresponding to different stages of the analysis.

#### Usage of `pipeline_aDNA.py`

The `pipeline_aDNA.py` file is the main entry point for the aDNA pipeline. It orchestrates the execution of various scripts and tasks, guiding the pipeline through its different stages.

##### Pipeline Stages

The pipeline is divided into several stages, executed sequentially:

1. **Raw reads processing** - for more details see [Raw Read Processing](raw_reads_processing.md)
2. **Reference genome processing** - for more details see [Reference Genome Processing](ref_genome_processing.md)
3. **Additional analysis**   

The pipeline automatically manages dependencies and workflow execution.

##### Running the Pipeline

To run the pipeline, navigate to the root directory containing the scripts folder and execute:

```bash
python scripts/pipeline_aDNA.py
```

#### Notes for running the Pipeline

##### Running the Pipeline in the Background

Depending on the size of the data, it may take some time to complete the pipeline. Thus it is recommended to run the pipeline in the background. You can do this by running the following command:

```bash
nohup python -u scripts/pipeline_aDNA.py > pipeline.log 2>&1 &
```

##### Restarting the Pipeline

The pipeline will always start from the first stage, even if it was previously completed. The individual steps recognize the state of the pipeline and will start from the last completed stage. Completed steps will be skipped. 

If you want to restart the pipeline from the beginning, you can delete the relevant folders and re-run the pipeline. This will lead to a complete re-processing of the data.

##### Parallelization

Some stages support parallelization. The number of threads can be adjusted in the config file file.

## Requirements

There are various dependencies which are required to run the aDNA pipeline. They can be installed as a conda environment.

### Step 1: Get environment file

* [adna_pipeline_env.yml](adna_pipeline_env.yml)

This file contains all bio informatic tools (fastp, angsd, ..), python libraries (pandas, pysam, ...) and R libraries (ggplot2, dplyr, ... ) for the pipeline itself as well as for all the used programs and tools.

### Step 2: Create the environment

```bash
conda env create -f adna_pipeline_env.yml
```

### Step 3: Activate it

```bash
conda activate adna_pipeline
```
