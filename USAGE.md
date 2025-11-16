# aDNA Pipeline — Walkthrough & Usage

This document walks through how the `aDNA_Pipeline_Snakemake` workflow works and how to run it. It complements the repository `README.md` by giving a concise, practical guide to getting the pipeline running, configuring it, and understanding the main stages and outputs.

## Table of contents

- Quick summary
- Prerequisites
- Where to run
- Project layout (relevant paths)
- Input file naming & species folders
- How the pipeline is structured (high level)
- Configuring the pipeline
- Running the pipeline — common commands
- Conda environments and external tools
- Outputs and where to find results
- Troubleshooting & common problems
- Examples
- Development & adding rules
- Where to look for more detail


## Quick summary

- Workflow engine: Snakemake (minimum version: 9.9.0 as enforced by the pipeline).
- Main entrypoint: `workflow/Snakefile` (a thin wrapper that includes rule files).
- Config file: `config/config.yaml` controls pipeline behavior, stages, species, and tool-specific settings.
- Major stages: raw reads processing, reads → reference processing (mapping, damage), and additional analyses (contamination checks: ECMSD, Kraken, Centrifuge).

## Prerequisites

- Linux or macOS (this workspace was developed on macOS/Unix-like systems).
- Snakemake (>= 9.9.0). Install via conda or pip:

```bash
# Recommended: create a conda env and install snakemake
conda create -n snakemake python=3.11 snakemake
conda activate snakemake
```

- Conda for `--use-conda` mode (the pipeline uses rule-level conda environment YAMLs).
- Tool-specific binaries/databases when using contamination tools (centrifuge, kraken, ecmsd). The paths and conda env YAMLs are referenced in `config/config.yaml` and each rule's settings.

## Where to run

The Snakemake `Snakefile` lives in `workflow/Snakefile`. You can either:

- cd into the project root (where `workflow/` is located) and run snakemake (Snakemake will load `workflow/Snakefile` because the top-level Snakefile imports it). Or
- call snakemake explicitly with `-s workflow/Snakefile` from the project root.

Both approaches are valid. Examples below assume you're in the project root (`aDNA_Pipeline_Snakemake/`).

## Project layout (relevant paths)

- `workflow/Snakefile` - main pipeline entry point (includes files under `workflow/rules/` and `workflow/scripts/`).
- `config/config.yaml` - pipeline config with global/pipeline/species settings.
- `workflow/rules/` - rule collections for each processing stage (raw reads, mapping, analytics, plotting).
- `workflow/scripts/` - helper Python/Snakemake scripts used by rules.
- `reports/` - generated reports and workflow provenance (`reports/workflow.rst` is referenced by the pipeline).
- `<species>/raw/reads/` - place raw FASTQ.gz reads here (see filenames below).
- `<species>/raw/ref/` - reference files for mapping.
- `<species>/processed/` and `<species>/results/` - created by the pipeline for intermediate and final outputs.

## Input file naming & species folders

The pipeline expects inputs to be placed under a species folder with the following structure (examples):

```
<species>/
  raw/
    reads/            # raw fastq.gz files
    ref/              # reference FASTA, indexes are created in processing
    mtdna/            # optional mitochondrial reads
```

Raw read filename convention expected by the pipeline:

```
<Individual>_<Original_Filename>.fastq.gz
```

Example: `Bger1_326862_S37_R1_001.fastq.gz`.

If you add a new species, add an entry under `species:` in `config/config.yaml` (the `species` keys map to folder names)

## How the pipeline is structured (high level)

The pipeline is split into two main grouped workflows which are included by `workflow/Snakefile`:

1. `raw_reads_processing.smk` (rules under `workflow/rules/raw_read/...`)
   - prepare raw reads
   - adapter removal
   - quality filtering
   - merging reads by individual
   - FastQC & MultiQC quality reports
   - contamination checks (ECMSD, Centrifuge, Kraken)
   - various read-level analytics & plotting

2. `reads_to_reference_processing.smk` (rules under `workflow/rules/reads_to_reference/...`)
   - prepare reference for mapping (indexing)
   - map reads to reference
   - analyze damage and rescale BAMs
   - coverage statistics & plotting
   - determine endogenous read counts

Control flags: each pipeline stage and many individual process steps have an `execute: true|false` flag in `config/config.yaml`. If a block is missing, the pipeline defaults to `execute: true`.

## Configuring the pipeline

Open `config/config.yaml`. Key sections:

- `project_name`: human-friendly name for reports.
- `pipeline:`: top-level stages (`raw_reads_processing`, `reference_processing`, etc.). Each stage has `execute:` and nested process steps with their own `execute` flags and `settings`.
- `species:`: one entry per species (folder name), with `name:` used for human-readable labels.
- Tool-specific settings for contamination tools include `conda_env`, `executable`, and `database` paths.

To disable a stage, set `pipeline.<stage>.execute: false`. To disable a specific process step, set that step's `execute: false`.

Example: to only run mapping and downstream analyses but skip raw-read processing, set in `config/config.yaml`:

```yaml
pipeline:
  raw_reads_processing:
    execute: false
  reference_processing:
    execute: true
```

## Running the pipeline — common commands

Dry-run (no execution, shows what would run):

```bash
snakemake --cores 4 --use-conda -n -s workflow/Snakefile
```

Full run (resume from previous state):

```bash
snakemake --cores 8 --use-conda --keep-going -s workflow/Snakefile
```

Run in background and capture logs:

```bash
nohup snakemake --cores 8 --use-conda --keep-going -s workflow/Snakefile > pipeline.log 2>&1 &
```

Run a single rule / target (for debugging):

```bash
# Example: only build the MultiQC report or another specific output
snakemake --cores 2 --use-conda path/to/target_file -s workflow/Snakefile
```

Run with increased verbosity and job provenance (helpful for debugging):

```bash
snakemake --cores 4 --use-conda --verbose --printshellcmds -s workflow/Snakefile
```

Tips for large clusters or cluster submission:
- Use Snakemake's cluster/DRMAA support or `--cluster`/`--kubernetes` flags.
- Provide `--jobs` / `--cluster` configurations and per-rule `resources` if needed.

## Conda environments and external tools

- Many rules declare rule-level conda envs (rule `conda:`). When using `--use-conda`, Snakemake will create and activate these envs automatically.
- `config/config.yaml` references some conda env YAMLs for contamination tools (example: `../../../../envs/centrifuge.yaml`). Verify these paths exist or replace them with local env YAMLs.
- Some tools (ECMSD, Kraken, Centrifuge) require external databases; configure their `database:` paths in `config/config.yaml`.

## Outputs and where to find results

- Intermediate files: `<species>/processed/` (per-step intermediates)
- Final results & reports: `<species>/results/`
- Plots and analytics: located under species results and `reports/`
- Workflow provenance (summary) is collected in `reports/workflow.rst` and rule-level logs are printed to stdout unless redirected.

## Troubleshooting & common problems

- ECMSD failures: ECMSD calls can fail for some inputs — the pipeline authors know this. Use `--keep-going` to continue with other rules when ECMSD fails, or disable ECMSD in `config/config.yaml`:

```yaml
pipeline:
  raw_reads_processing:
    contamination_analysis:
      tools:
        ecmsd:
          execute: false
```

- Missing conda env YAMLs: check rule declarations under `workflow/rules/*` or `workflow/scripts/` and ensure the referenced env YAMLs exist and are readable.
- Missing databases: contamination steps require local database files. Update their `database` paths in `config/config.yaml`.
- Permissions / disk space: mapping, indexing, and intermediate files can use a lot of disk. Ensure enough space and correct write permissions for species folders.

If a rule fails, Snakemake prints the failing shell command. Re-run with `-n` to inspect, or `--rerun-incomplete` if outputs are partially written and you want Snakemake to re-run them.

## Quick checklist before running

- Confirm species folders exist and raw reads are placed in `<species>/raw/reads/` using the filename convention.
- Update `config/config.yaml` with correct `project_name`, `species` entries, and any tool/database paths.
- Verify conda environments referenced in `config.yaml` or rule `conda:` statements are available (or will be created by Snakemake).
- Ensure sufficient disk space and write permissions for `processed/` and `results/` folders.
- If using heavy tools (centrifuge/kraken), confirm database paths are correct and reachable.

## Common Snakemake commands (cheat sheet)

- Dry-run (shows what would run):

```bash
snakemake --cores 4 --use-conda -n -s workflow/Snakefile
```

- Run full pipeline with 8 cores and keep going on errors:

```bash
snakemake --cores 8 --use-conda --keep-going -s workflow/Snakefile
```

- Run a single target (useful for debugging):

```bash
snakemake --cores 2 --use-conda path/to/target_file -s workflow/Snakefile
```

- Re-run incomplete/partial outputs (force re-execution of incomplete outputs):

```bash
snakemake --cores 4 --use-conda --rerun-incomplete -s workflow/Snakefile
```

- Use a submission cluster (example pattern using GNU parallel / local cluster):

```bash
snakemake --jobs 100 --cluster "qsub -V -cwd -l h_vmem=8G" -s workflow/Snakefile
```

## Tips & troubleshooting commands

- Increase per-rule threads/resources: edit `config/config.yaml` or adjust rule `threads:`/`resources:` in the rule files under `workflow/rules/` if needed. Many rules read default thread counts from `config`.
- Disable noisy or failing steps (quick): set the specific `execute: false` flag in `config/config.yaml` for the offending stage (for example, disable ECMSD if it fails frequently).
- Inspect the last failing command and stderr/stdout: run Snakemake without `--quiet` and look at the printed shell command; rule logs (if created) are stored next to outputs or printed to stdout.
- Show detailed traceback for a rule failure:

```bash
snakemake --cores 1 --use-conda --keep-going -s workflow/Snakefile --printshellcmds --verbose
```

- Re-create a conda env manually if automatic creation fails: find the rule's `conda:` YAML (search `workflow/rules/` for `conda:`) and run `conda env create -f path/to/env.yaml -n some_env_name`, then set `CONDA_DEFAULT_ENV` or run Snakemake with `--conda-prefix` pointing to an envs directory.

## Small maintenance & developer notes

- To add a new species, add an entry under `species:` in `config/config.yaml` and create the folder structure under the repository root (or in your project workspace) with `raw/reads/` and `raw/ref/`.
- When adding new rules keep them organized in `workflow/rules/` and include them in the appropriate group file (`raw_reads_processing.smk` or `reads_to_reference_processing.smk`).
- For reproducibility, commit `config/config.yaml` changes to a branch, but keep machine-specific absolute database paths out of shared config or document them clearly.


## Examples

1) Run a quick dry-run for 4 threads:

```bash
snakemake --cores 4 --use-conda -n -s workflow/Snakefile
```

2) Run full pipeline, keep going on errors, 8 threads:

```bash
snakemake --cores 8 --use-conda --keep-going -s workflow/Snakefile
```

## Development & adding rules

- The pipeline is modular. Rule files live under `workflow/rules/` and are included by grouped files like `raw_reads_processing.smk` and `reads_to_reference_processing.smk`.
- New rules should follow the repository's pattern: place them in a logical subfolder and include them from the appropriate group file.
- Use `config/config.yaml` to add or toggle settings instead of hardcoding values in rules.

## Where to look for more detail

- `README.md` (top-level) — conceptual overview (already present in this repo).
- `config/config.yaml` — pipeline settings and species list.
- `workflow/rules/` and `workflow/scripts/` — actual execution logic and helper scripts.

## Contact / Contributions

If something is unclear or you want help adapting the pipeline to your data, open an issue or contact the repository owner.

---

Last updated: 22 October 2025
