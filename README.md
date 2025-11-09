# aDNA Pipeline

This project contains a pipeline to analyze raw ancient data, obtained from the sequencing facility. The pipeline includes various Snakemake workflows to process, analyze, and generate reports on the sequence quality, which helps decide if an aDNA extraction and sequencing was successful, and further polishes the data for downstream analyses.

Note: This pipeline is still in development, but can already be used for analysis.

## Setup Overview

- Before running the pipeline, ensure you have an environment with Snakemake and the required dependencies installed.
- Required dependencies for pipeline processing will be installed automatically, except for the contamination analysis tools. Those have to be installed separately and their details need to be added to the config file.
- You need to add the project path to the config file.
- you need to add species details to the pipeline.

When adding a new species, make sure to 
- the species folder should be placed in the root folder of your pipeline
- add the folder name to the `config.yaml` below `species:` 
- your reads should be renamed according to the naming convention specified in the [RAW Reads filenames](#RAW-Reads-filenames) section
- the pipeline supports automatically moving the raw reads to the `<species>/raw/reads/` folder as well as the reference genome to the `<species>/raw/ref_genome/` folder. Simply provide the files in the `<species>` folder. Alternatively, you can manually move the files to the respective folders.
  - provide the raw reads in `<species>/raw/reads/` folder
  - provide the reference genome in `<species>/raw/ref_genome/` folder
- all other folders will be created and populated automatically
  - folder `<species>/processed/` contains the intermediary files during processing. Most of these files are marked as temporary and will be deleted at the end of the pipeline. Some files are kept to allow reprocessing the pipeline from different points in case something fails.
  - folder `<species>/results/` contains the final results and reports. 
  - general reads processing data will be in either `processed`or `results`. Everything related to a reference genome will have a `<reference_genome>` folder under `processed`or `results`. Typically, only the `results` folder will contain information required for further analyis. In case more information is required, the original files can often be found in the `processed` folder. Some exemptions include `*.sam` and unsorted `*.bam` files. These are deleted to save storrage space. Most other files are kept in order to allow reprocessing the pipeline from different points in case something fails. If a step should be repeated, the relevant files need to be deleted manually. 

## Configuration File Structure for aDNA Pipeline (`config.yaml`)

The `config.yaml` file is used to configure the aDNA pipeline. It contains settings such as project name, the species list and the pipeline stages and process steps.

### Global Settings

* **project\_name**: Name of the aDNA project.

### Pipeline Settings

Defines the overall pipeline behavior, including execution controls and process details.

#### Pipeline Stages and Process Steps

* The pipeline is broken into **stages** (e.g., `raw_reads_processing`, `reference_genome_processing`, `post_processing`).
* Each stage contains multiple **process steps** (e.g., `adapter_removal`, `deduplication`, ...).
* Both stages and process steps can be controlled with `execute: true/false` flags to enable or disable them.
* Some process steps include additional configurable settings (e.g., adapter sequences, database paths, ...).

#### Important Defaults

* You **do not need to specify all stages or process steps** explicitly.
* Any **stage or process step not provided in the config defaults to `execute: true`** and will be executed.

### Example `config.yaml`

```yaml
# config.yaml - Configuration file for aDNA pipeline

# Global settings
project_name: "aDNA_Project"

pipeline:
  raw_reads_processing:     # stage
    execute: true           # disable or enable this stage
    quality_checking_raw:   # process step
      execute: true         # disable or enable this process step
        adapter_removal:
      execute: true
      settings: 
        min_quality: 5    # default 5
        min_length: 15    # default 15
        adapters_sequences: # if not provided, fastp will try to identify adapters
          r1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
          r2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" 
    quality_checking_trimmed:
      execute: true
    quality_filtering:
      execute: true
      settings:
        min_quality: 15   # default 15
        min_length: 15    # default 15
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
            conda_env: "../../../../envs/centrifuge.yaml"
            executable: "centrifuge"
            database: "/mnt/data5/sarah/aDNA/centrifuge_db"
        kraken: 
          execute: true
          settings:
            conda_env: "../../../../envs/kraken.yaml"
            executable: "kraken"
            database: "/mnt/data5/sarah/aDNA/kraken_db"
        ecmsd:
          execute: true
          settings:
            conda_env: "../../../../envs/ecmsd.yaml"
            executable: "/mnt/data2/sarah/app_ecmsd/shell/ECMSD.sh"
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
* The `--keep-going` flag allows the pipeline to continue even if a rule fails. This is useful for debugging purposes. Somtimes the analysis of ECMSD fails due to issues with the input data. In this case, the rest of the pipeline can still be executed.
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

### Output Folder Structure

```
<species>/
├── raw/
│   ├── ref_genome/
│   └── reads/
├── results/
│   ├── <ref_genome>/
│   │   ├── endogenous/
│   │   │   ├── <Individual>/
│   │   ├── damage/
│   │   │   └── <Individual>/
│   │   ├── statistics/
│   │   │   └── <Individual>/
│   │   ├── coverage/
│   │   │   └── <Individual>/
│   │   └── plots/
│   │       ├── endogenous_reads/
│   │       └── coverage/
│   ├── contamination_analysis/
│   │   └── ecmsd/
│   │       └── <read_name>/
│   └── reads/
│       ├── reads_merged/
│       │   └── fastqc/
│       ├── reads_raw/
│       │   └── fastqc/
│       ├── reads_quality_filtered/
│       │   ├── fastqc/
│       │   └── fastp_report/
│       ├── reads_trimmed/
│       │   └── fastqc/
│       └── statistics/
├── logs/
│   ├── <ref_genome>/
│   │   ├── endogenous/
│   │   │   └── <Individual>/
│   │   ├── consensus/
│   │   │   └── <Individual>/
│   │   ├── damage/
│   │   │   └── <Individual>/
│   │   ├── statistics/
│   │   │   └── <Individual>/
│   │   ├── coverage/
│   │   │   └── <Individual>/
│   │   └── plots/
│   │       ├── endogenous_reads/
│   │       └── coverage/
│   ├── contamination_analysis/
│   │   └── ecmsd/
│   └── reads/
│       ├── reads_merged/
│       │   └── fastqc/
│       ├── reads_raw/
│       │   └── fastqc/
│       ├── reads_quality_filtered/
│       │   └── fastqc/
│       └── reads_trimmed/
│           └── fastqc/
├── processed/
│   ├── <ref_genome>/
│   │   ├── consensus/
│   │   │   └── <Individual>/
│   │   └── mapped/
│   └── reads/
│       └── reads_merged/
```