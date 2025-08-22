# Processing reads and mapping them to a reference genome

This pipeline processes ancient DNA (aDNA) reads by mapping them to a reference genome, analyzing damage, calculating coverage metrics, and producing visual summaries. Below is a step-by-step breakdown of the process:

## Process Overview

### 1. **Prepare Reads for Mapping**

* **Per Individual and Per Species FASTQ Creation**
  Reads from each individual are processed to generate one FASTQ.GZ file per individual.

### 2. **Prepare Reference Genome**

* Reference genomes are preprocessed for mapping (e.g., indexing).

### 3. **Map Reads to Reference Genome**

* Processed FASTQ files are aligned to the corresponding reference genome.

> 📝 *Note (2025-04-07): SAM to BAM conversion, sorting, and indexing now occur immediately after mapping to save disk space. This step is no longer run separately.*

### 4. **Analyze DNA Damage**

* `mapDamage` is run on mapped reads to assess damage patterns typical for aDNA (e.g., deamination at read ends).

### 5. **Quality Control: Endogenous Content**

* Endogenous reads are quantified per individual.

### 6. **Determine Coverage Depth and Breadth**

* **Coverage metrics are calculated per individual**:

  * **Average depth**
  * **Maximum depth**
  * **Covered vs. total bases**
  * **Percent covered**
* Results are aggregated into summary and detailed reports for further analysis and plotting.

### 7. **Generate Summary Plots**

* Visualizations are produced for:

  * **Coverage depth and breadth**
  * **Endogenous DNA content**
* Plots are generated per individual and per species to support comparative analysis.

> 🛑 *Optional step (commented out in code):*
> Special sequence extraction (`extract_special_sequences`) can be included if needed for downstream analysis.

## Additional Info: Analysis

### Coverage Depth and Breadth Analysis

This analysis quantifies the depth and breadth of coverage of ancient DNA (aDNA) sequencing data mapped to a reference genome. The analysis is performed at two levels: **individual** and **combined (species-level)**.

#### Data Structure

##### Individual Coverage Outputs
Found in:
`aDNA/<species>/processed/<ref_genome>/coverage_depth_breadth/<individual>`

Each individual has a `.tsv` output from `samtools depth` which is used as a base for the summary files below.

##### Individual Summary CSV
Found in:
`aDNA/<species>/results/<ref_genome>/coverage_depth_breadth`
File: `<individual>_<ref_genome>_extended_coverage_analysis.csv`

This file contains **scaffold-level coverage statistics** for a single individual. It is used to summarize mapping quality across the reference genome and is suitable for creating **summary barplots** per individual.

###### File Structure

Each row represents a **single scaffold** from the reference genome, with calculated coverage metrics based on the individual's mapped reads.

**Columns:**

| Column Name       | Description                                                                        |
| ----------------- | ---------------------------------------------------------------------------------- |
| `scaffold`        | Name/ID of the scaffold or contig in the reference genome.                         |
| `avg_depth`       | Average sequencing depth across the scaffold.                                      |
| `max_depth`       | Maximum depth observed at any base in the scaffold.                                |
| `covered_bases`   | Number of scaffold bases with at least one read mapped.                            |
| `total_bases`     | Total number of bases in the scaffold.                                             |
| `percent_covered` | Proportion of the scaffold covered by reads (`covered_bases / total_bases * 100`). |

> 🔹 **Note**: Unlike the detailed file, this version contains no `Filename` column since the file itself is specific to one individual.

###### Example Entry

```
scaffold,avg_depth,max_depth,covered_bases,total_bases,percent_covered
JPZV02000014.1,1.0854,137,840600,1402747,59.93
```

This line indicates that:

* Scaffold `JPZV02000014.1` has an average depth of \~1.09× and a maximum depth of 137.
* \~60% of its bases are covered by reads from this individual (`Bger1`).

##### Detailed Combined CSV (Per Scaffold and Individual)**
Found in:
`aDNA/<species>/results/<ref_genome>/coverage_depth_breadth`
File: `<species>_combined_coverage_analysis_detailed.csv`

This file contains **scaffold-level coverage statistics per individual**. It is used for detailed visualizations such as **violin plots** to assess coverage variation across scaffolds and samples.

###### File Structure

Each row represents a single **scaffold** for a given individual.

**Columns:**

| Column Name       | Description                                                                                    |
| ----------------- | ---------------------------------------------------------------------------------------------- |
| `scaffold`        | Name/ID of the scaffold or contig in the reference genome.                                     |
| `avg_depth`       | Average sequencing depth across the scaffold (mean number of reads covering each base).        |
| `max_depth`       | Maximum depth observed at any single base on the scaffold.                                     |
| `covered_bases`   | Number of bases in the scaffold with at least one read mapped (coverage > 0).                  |
| `total_bases`     | Total number of bases in the scaffold.                                                         |
| `percent_covered` | Percentage of the scaffold covered by at least one read (`covered_bases / total_bases * 100`). |
| `Filename`        | The name of the BAM file corresponding to the individual, identifying the sample.              |

###### Example Entry

```
scaffold,avg_depth,max_depth,covered_bases,total_bases,percent_covered,Filename
KZ614359.1,0.2674,281,218910,1214666,18.02,Bger2.fastq_GCA_000762945.2_Bger_2.0_genomic_sorted.bam
```

This example shows that:

* Scaffold `KZ614359.1` has low average depth (`0.27×`) with a max depth of 281.
* Only \~18% of the scaffold bases are covered by sequencing reads.
* The data corresponds to the individual sample `Bger2`.

##### Combined Summary CSV (Species Level)
Found in:
`aDNA/<species>/results/<ref_genome>/coverage_depth_breadth`
File: `species_combined_coverage_analysis.csv`

This file provides a **summary of coverage statistics per individual**, aggregated across all scaffolds in the reference genome. It is useful for comparing overall sequencing quality between individuals and is typically used in **barplots** or summary tables.

###### File Structure

Each row represents a single **individual**, identified by their BAM filename, with global coverage metrics calculated over the entire reference genome.

**Columns:**

| Column Name             | Description                                                      |
| ----------------------- | ---------------------------------------------------------------- |
| `Filename`              | BAM file name corresponding to the individual sample.            |
| `OverallAvgDepth`       | Average depth across the entire reference genome.                |
| `OverallMaxDepth`       | Maximum read depth observed at any position in the genome.       |
| `OverallCoveredBases`   | Total number of bases with at least one mapped read.             |
| `OverallTotalBases`     | Total number of bases in the reference genome.                   |
| `OverallPercentCovered` | Percentage of the reference genome covered by at least one read. |

###### Example Entry

```
Filename,OverallAvgDepth,OverallMaxDepth,OverallCoveredBases,OverallTotalBases,OverallPercentCovered
Bger2.fastq_GCA_000762945.2_Bger_2.0_genomic_sorted.bam,21.01,3033736,1672203996,1791714300,93.33
```

This row shows that:

* The individual `Bger2` achieved an average genome-wide depth of \~21.
* \~93% of the reference genome was covered by at least one read.
* The maximum depth at any position was over 3 million reads.

#### Analysis Summary

For each **individual**:

* Coverage is calculated **per scaffold**.
* Metrics include average depth, max depth, total and covered bases, and percent of scaffold covered.
* Summary metrics are generated by computing total or average values across all scaffolds.

For **combined analyses**:

* Coverage values are merged across individuals.
* Suitable for comparisons across samples or species-level overviews.

#### Visualization

* **Barplots**: Use individual summary data (`<species>_coverage_analysis.csv`) to compare total coverage.
* **Violin plots**: Use detailed data (`<species>_coverage_analysis_detailed.csv`) to show distribution of coverage across scaffolds and individuals.