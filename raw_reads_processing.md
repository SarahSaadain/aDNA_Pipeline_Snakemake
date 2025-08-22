# Raw Read Processing

This pipeline part processes raw ancient DNA (aDNA) sequencing reads through quality control, adapter removal, filtering, deduplication, and contamination analysis. It produces both summary statistics and diagnostic plots to assess data quality before downstream analyses.

## Step-by-Step Process


### 1. **Initial Quality Control**

* **Tools**: `FastQC`, `MultiQC`
* Perform QC on raw reads to assess quality scores, GC content, adapter contamination, and sequence length distribution.

### 2. **Adapter Removal and Read Merging**

* **Tool**: `fastp`
* Remove sequencing adapters from both single-end and paired-end reads.
* Merge overlapping paired-end reads to improve downstream alignment.

Script: `scripts/raw_reads_processing/execute_fastp_adapter_remove_and_merge.py`
* Parameters:
	+ Paired-end reads:
		- `-i`: Input file (R1 and R2)
		- `-o`: Output file (merged)
		- `--adapter_sequence`: Adapter sequence
		- `--length_required`: Minimum length of reads to keep (set to `15`)
		- `--trim_poly_x`: Trim poly-X tails (set to `5`)
		- `--qualified_quality_phred`: Minimum quality score to keep (set to `5`)
		- `--unqualified_percent_limit`: Maximum percentage of unqualified bases (set to `40`)
		- `--n_base_limit`: Maximum number of N bases allowed (set to `5`)
	+ Single-end reads:
		- `-i`: Input file
		- `-o`: Output file
		- `--adapter_sequence`: Adapter sequence (not specified)
		- `--length_required`: Minimum length of reads to keep (set to `15`)
		- `--trim_poly_x`: Trim poly-X tails (set to `5`)
		- `--qualified_quality_phred`: Minimum quality score to keep (set to `5`)
		- `--unqualified_percent_limit`: Maximum percentage of unqualified bases (set to `40`)
		- `--n_base_limit`: Maximum number of N bases allowed (set to `5`)

### 3. **QC After Adapter Removal**

* **Tools**: `FastQC`, `MultiQC`
* Validate adapter removal success and check for improvements in read quality.

### 4. **Quality Filtering**

* **Tool**: `fastp`
* Apply base-level quality filtering to remove low-quality reads or read segments.

Script: `scripts/raw_reads_processing/execute_fastp_quality_filter.py`

* Parameters:
	+ `-i`: Input file
	+ `-o`: Output file
	+ `-q`: Minimum quality score to keep (set to `30`)
	+ `-p`: Minimum percent of bases that must have `-q` quality (set to `75`)
	+ `--length_required`: Minimum length of reads to keep (set to `15`)
	+ `--trim_poly_x`: Trim poly-X tails (set to `5`)
	+ `--qualified_quality_phred`: Minimum quality score to keep (set to `5`)
	+ `--unqualified_percent_limit`: Maximum percentage of unqualified bases (set to `40`)
	+ `--n_base_limit`: Maximum number of N bases allowed (set to `5`)

### 5. **QC After Quality Filtering**

* **Tools**: `FastQC`, `MultiQC`
* Evaluate the effectiveness of quality filtering.

### 6. **Duplicate Removal**

* **Tool**: `fastp`
* Remove PCR duplicates to retain unique reads and reduce artificial coverage inflation.

Script: `scripts/raw_reads_processing/execute_fastp_deduplication.py`

* Parameters:
	+ `-i`: Input file
	+ `-o`: Output file
	+ `--dedup`: Enable deduplication (set to `True`)
	+ `--length_required`: Minimum length of reads to keep (set to `15`)
	+ `--trim_poly_x`: Trim poly-X tails (set to `5`)
	+ `--qualified_quality_phred`: Minimum quality score to keep (set to `5`)
	+ `--unqualified_percent_limit`: Maximum percentage of unqualified bases (set to `40`)
	+ `--n_base_limit`: Maximum number of N bases allowed (set to `5`)

### 7. **QC After Deduplication**

* **Tools**: `FastQC`, `MultiQC`
* Assess final read quality post-deduplication.

### 8. **Quality Report Generation**

* An HTML report is generated to summarize all quality control results from each step of the pipeline.

### 9. **Read Processing Statistics**

* Summarize read counts before and after each processing step.
* Useful for evaluating data retention and filtering efficiency.

### 10. **Read Length Distribution**

* Analyze the distribution of read lengths to identify fragmentation patterns typical of ancient DNA.

### 11. **Contamination Detection**

* **Tools**: `Centrifuge` and `Kraken`
* Classify and quantify possible contaminant sequences (e.g., microbial or modern DNA).

### 12. **Summary Plot Generation**

* Visualize:

  * Read counts before/after each processing step
  * Sequence length distributions plots are generated for raw reads, adapter-removed reads, quality-filtered reads, and deduplicated reads.
  * Contamination profiles

This pipeline ensures that raw aDNA data is cleaned, filtered, and verified for quality before alignment to reference genomes or further analysis.

## FastQC and MultiQC

FastQC and MultiQC are used to generate quality control reports for raw reads, adapter-removed reads, quality-filtered reads, and deduplicated reads. The reports are saved in the `species/results/quality_control/fastqc` and `species/results/quality_control/multiqc` folders, respectively.
