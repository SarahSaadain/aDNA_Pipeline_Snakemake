# TEplotter

A Python toolkit for converting BAM/SAM alignment files to sequence overview (SO) format, with support for variant detection (SNPs, indels), coverage normalization, and visualization-ready outputs.

## Overview

TEplotter processes short-read alignments to extract coverage, SNP, and indel information for each reference sequence. It's particularly useful for analyzing transposable element (TE) coverage and variation across genomic regions.

### Workflow

```
BAM/SAM → bam2so.py → SO-file → normalize-so.py → So-file 
                      ↓                              ↓
                       →                              → so2plotable.py → plotable-file → visualize-plotable.R → png,eps,pdf,svg


# estimate copy numbers of TEs or any other sequence of interest
BAM/SAM → bam2so.py → SO-file → normalize-so.py → SO-file -> estimate-SO.py -> copy-number-estimates
```

## Features

- **BAM/SAM Processing**: Converts alignments to sequence overview (SO) format with pysam
- **Variant Detection**: Identifies SNPs and indels with configurable thresholds
- **Coverage Analysis**: Tracks read depth per position with ambiguous mapping support
- **Normalization**: Normalizes coverage using single-copy genes (SCGs) as reference
- **Copy Number Estimation**: Computes average copy number from coverage data
- **Visualization**: Generates tab-delimited output compatible with R/ggplot2

## Installation

### Requirements

- Python 3.7+
- pysam
- (Optional) samtools (for FASTA indexing)

### Quick Install

```bash
pip install pysam
```

## Usage

### 1. Convert BAM/SAM to Sequence Overview Format

```bash
python bam2so.py \
  --infile alignments.bam \
  --fasta reference.fasta \
  --mapqth 5 \
  --mc-snp 5 \
  --mf-snp 0.1 \
  --mc-indel 3 \
  --mf-indel 0.01 \
  --output-file output.so
```

**Parameters:**
- `--infile`: Input BAM or SAM file (required)
- `--fasta`: Reference FASTA file (required)
- `--mapqth`: Mapping quality threshold (default: 5) — reads below this are "ambiguous"
- `--mc-snp`: Minimum SNP count (default: 5)
- `--mf-snp`: Minimum SNP frequency (default: 0.1)
- `--mc-indel`: Minimum indel count (default: 3)
- `--mf-indel`: Minimum indel frequency (default: 0.01)
- `--output-file`: Output SO file; if omitted, prints to stdout

### 2. Normalize Coverage by Single-Copy Genes

```bash
python normalize-so.py \
  --so input.so \
  --scg-end _scg \
  --end-distance 100 \
  --exclude-quantile 25 \
  --output-file normalized.so
```

**Parameters:**
- `--so`: Input SO file (required)
- `--scg-end`: Suffix to identify single-copy genes (default: "_scg")
- `--end-distance`: Distance from sequence ends to exclude from normalization (default: 100)
- `--exclude-quantile`: Quantile threshold for excluding extreme coverage (default: 25)
- `--output-file`: Output file; if omitted, prints to stdout

### 3. Estimate Copy Number

```bash
python estimate-so.py \
  --so normalized.so \
  --end-distance 100 \
  --exclude-quantile 25 \
  --output-file copy_numbers.txt
```

Outputs tab-delimited format: `seqname <tab> average_coverage <tab> length`

### 4. Convert to Plotable Format

```bash
python so2plotable.py \
  --so input.so \
  --sample-id year1933 \
  --seq-ids gypsy,act_scg \
  --output-file sequences.plotable or --output-dir myplotables
```

Generates visualization-ready tab-delimited output with columns:
- seqname, sampleid, feature (cov/ambcov/snp/del/ins), position, value
  
**Note:** normalized and not-normalized SO files may be visualized! The idea is that the unmodified raw-data may be used as well as the normalized coverage adjusted to the number of reads 
**Note:** its possible to provide '--seq-ids all', in which case all sequences will be prepared for plotting (i.e. converted into the plotable format)
**Note:** an important design decision was that plotable entries for different samples (e.g strains collected at different years) and sequences (e.g. TEs or SCGs) can be combined freely by the user, which allows for joint visualization

### 5. Visualize

```bash
Rscript visualize-plotable.R sequences.plotable output.png
```
Any other extension may be used. Importantly plotable-files from different samples may be concatenated (using cat) which will automatically invoke facetting. this allows to test in which samples copy number of TEs change

## File Formats

### Sequence Overview (SO) Format

A tab-delimited format containing:
- Coverage per position
- Ambiguous coverage (low MAPQ)
- SNP calls with allele counts
- Indel calls (insertions/deletions) with counts

### Plotable Format

Tab-delimited format for direct plotting in R/ggplot2:
```
seqname  sampleid  feature  position  value
TE_001   sample_1  cov      1         42.5
TE_001   sample_1  snp      100       A→T:3
```

## Examples

### Complete Pipeline

```bash
# 1. Convert BAM to SO
python bam2so.py \
  --infile reads.bam \
  --fasta te_library.fasta \
  --output-file raw.so

# 2. Normalize coverage
python normalize-so.py \
  --so raw.so \
  --scg-end _scg \
  --output-file normalized.so

# 3. Estimate copy numbers
python estimate-so.py \
  --so normalized.so \
  --output-file copy_numbers.txt

# 4. Format for plotting
python so2plotable.py \
  --so normalized.so \
  --sampleid my_sample \
  --output-file plotable.txt

# 5. Plot
Rscript visualize-plotable \
   plotable.txt \
   output.png
```

### Adjusting Variant Thresholds

Stricter SNP filtering:
```bash
python bam2so.py \
  --infile reads.bam \
  --fasta te_library.fasta \
  --mc-snp 10 \
  --mf-snp 0.2 \
  --output-file output.so
```

## Core Modules

### `modules.py`

Contains shared utilities:
- **SeqEntry**: Data class for sequence overview records
- **SeqBuilder**: Accumulates reads for a sequence and extracts variants
- **SeqEntryReader**: Iterator over SO files (supports gzip)
- **Writer**: Output writer (file or stdout)
- **NormFactor**: Normalization using SCG coverage
- **SNP, Indel**: Variant data classes

### `bam2so.py`

Main conversion pipeline:
1. Loads reference FASTA
2. Iterates alignments (skips unmapped, secondary, supplementary)
3. Builds per-position coverage, SNPs, indels
4. Writes SO format output

### `normalize-so.py`

Normalizes coverage using median coverage of SCGs across middle regions (excluding ends).

### `estimate-so.py`

Computes per-sequence average coverage (useful for copy number if normalized).

### `so2plotable.py`

Transforms SO to tab-delimited format suitable for visualization.

## Tips & Best Practices

1. **FASTA indexing**: Create a `.fai` index for faster access:
   ```bash
   samtools faidx reference.fasta
   ```

2. **BAM indexing**: Ensure BAM files are sorted and indexed:
   ```bash
   samtools sort -o sorted.bam input.bam
   samtools index sorted.bam
   ```

3. **Threshold tuning**: Adjust `--mc-snp`, `--mf-snp`, etc. based on coverage depth and expected variation rates.

4. **Single-copy genes**: If normalizing, ensure your reference includes sequences flagged with `_scg` suffix; coverage is median-normalized to these.

5. **Large files**: SO format can be space-intensive. Consider gzip compression for archival:
   ```bash
   gzip output.so
   # Later:
   python so2plotable.py --so output.so.gz ...
   ```

## Troubleshooting

### "Reference 'X' not found in FASTA"
Ensure the BAM header references sequences present in your FASTA file. Check consistency:
```bash
samtools view -H alignments.bam | head
samtools faidx reference.fasta | head
```

### Missing FASTA index
Create an index for faster reference lookups:
```bash
samtools faidx reference.fasta
```

### pysam import errors
Install or upgrade pysam:
```bash
pip install --upgrade pysam
```

## Author

Original author: Robert Kofler

## License

TODO

## Citation

If you use TEplotter in your research, please cite:
```
TODO
```

## Contributing

Contributions welcome! Please open an issue or submit a pull request.
