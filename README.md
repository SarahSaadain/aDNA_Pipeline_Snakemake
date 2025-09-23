# aDNA Pipeline

This project contains a pipeline to analyze raw ancient data, obtained from the sequencing facility. The pipeline includes various Snakemake workflows to process, analyze, and generate reports on the sequence quality, which helps decide if an aDNA extraction and sequencing was successful, and further polishes the data for downstream analyses.

Note: This pipeline is still in development, but can already be used for analysis.

## Setup Overview

- Before running the pipeline, ensure you have the necessary dependencies installed. Please refer to the [Requirements](#Requirements) section for the necessary dependencies and installation instructions.
- Raw reads and reference genomes must be provided in the relevant folders. 
    - Raw reads should be renamed according to the naming convention specified in the [RAW Reads filenames](#RAW-Reads-filenames) section. Also see the [Manually renaming the raw reads files](#Manually-renaming-the-raw-reads-files) section.
    - Reference genome must be provided in the `species/raw/ref_genome/` folder.
    - Species must be added to the config file.
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


### Example `config.yaml`

```yaml
# config.yaml - Configuration file for aDNA pipeline

# Global settings
project_name: "aDNA_Project"
project_description: "Analysis of ancient DNA data"

pipeline:
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
```

## Running the Pipeline

The aDNA pipeline is implemented using Snakemake, a workflow management system. Snakemake ensures reproducibility and efficient execution of the pipeline.

### Running the Pipeline

To run the pipeline, navigate to the root directory containing the `Snakefile` and execute:

```bash
snakemake --cores <number_of_threads> --use-conda
```

Replace `<number_of_threads>` with the number of CPU threads you want to allocate for the pipeline.

### Running the Pipeline in the Background

Depending on the size of the data, it may take some time to complete the pipeline. Thus it is recommended to run the pipeline in the background. You can do this by running the following command:

```bash
nohup snakemake --cores <number_of_threads> --use-conda > pipeline.log 2>&1 &
```

### Restarting the Pipeline

Snakemake automatically tracks the state of the pipeline and will only re-run steps that are incomplete or outdated. If you want to restart the pipeline from the beginning, you can delete the relevant output files and re-run the pipeline.

### Parallelization

The number of threads can be adjusted using the `--cores` option when running Snakemake. Additionally, the `threads` parameter in the `config.yaml` file can be used to set default thread usage for specific rules.

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


### Scripts Folder

The `scripts/` folder contains all necessary scripts for the aDNA pipeline, organized into subfolders corresponding to different stages of the analysis.

#### Usage of Snakemake Workflow

The `Snakefile` is the main entry point for the aDNA pipeline. It orchestrates the execution of various rules and tasks, guiding the pipeline through its different stages.

##### Pipeline Stages

The pipeline is divided into several stages, executed sequentially:

1. **Raw reads processing** - for more details see [Raw Read Processing](raw_reads_processing.md)
2. **Reference genome processing** - for more details see [Reference Genome Processing](ref_genome_processing.md)
3. **Additional analysis**

Snakemake automatically manages dependencies and workflow execution.

