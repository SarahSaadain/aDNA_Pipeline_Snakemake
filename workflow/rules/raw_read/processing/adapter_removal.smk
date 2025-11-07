def get_adapter_removal_input_reads(wc):
    """
    Returns a list of read files for R1/R2 if available.
    If only R1 exists, returns a single-element list [R1].
    """

    read_files = get_read_files_for_species(wc.species) 
    
    # turn read paths into file names only
    read_files = [os.path.basename(f) for f in read_files]

    reads_dir = f"{wc.species}/raw/reads"
 
    # R1
    base_r1 = f"{wc.sample}_R1"
    candidates_r1 = [f for f in read_files
                     if re.match(base_r1 + r"(\S*)?\.fastq\.gz", f)]
    if not candidates_r1:
        raise FileNotFoundError(f"No R1 found for {wc.sample}. Expected pattern: {base_r1}*.fastq.gz in {reads_dir}. Found files: {read_files}")

    r1 = os.path.join(reads_dir, sorted(candidates_r1)[0])
 
    # R2
    base_r2 = f"{wc.sample}_R2"
    candidates_r2 = [f for f in read_files
                     if re.match(base_r2 + r"(\S*)?\.fastq\.gz", f)]
 
    if candidates_r2:
        r2 = os.path.join(reads_dir, sorted(candidates_r2)[0])
        return [r1, r2]  # Paired-end
    else:
        return [r1]      # Single-end
 
# Rule: Adapter removal for single-end reads using fastp
rule fastp_se:
    input:
        sample=get_adapter_removal_input_reads,
    output:
        trimmed=temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.se.fastq.gz"),
        failed=temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.se.failed.fastq.gz"),
        html=report(
            "{species}/results/reads/reads_trimmed/fastp_report/{sample}_trimmed.se.html",
            #caption="../report/fastp.rst",
            category="quality control",
            subcategory="fastp",
        ),
        json="{species}/results/reads/reads_trimmed/fastp_report/{sample}_trimmed.se.json",
    message: "Trimming adapters from single-end reads in {input.sample}"
    log:
        "{species}/logs/reads/reads_trimmed/{sample}_trimmed.se.log",
    params:
        adapters=f"--adapter_sequence {config['pipeline']['raw_reads_processing']['adapter_removal']['settings']['adapters_sequences']['r1']}",
        extra=f"--length_required {config['pipeline']['raw_reads_processing']['adapter_removal']['settings'].get('min_length','15')} --trim_poly_x 5  --qualified_quality_phred {config['pipeline']['raw_reads_processing']['adapter_removal']['settings'].get('min_quality','5')} --unqualified_percent_limit 40 --n_base_limit 5"
    threads: workflow.cores
    wrapper:
        "v7.5.0/bio/fastp"
 
 
# Rule: Adapter removal for paired-end reads using fastp
rule fastp_pe:
    input:
        sample=get_adapter_removal_input_reads,
    output:
        trimmed=[
            temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.R1.fastq.gz"),
            temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.R2.fastq.gz"),
        ],
        # Unpaired reads separately
        unpaired1=temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.unpaired.R1.fastq.gz"),
        unpaired2=temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.unpaired.R2.fastq.gz"),
        merged=temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.fastq.gz"),
        failed=temp("{species}/processed/reads/reads_trimmed/{sample}_trimmed.pe.failed.fastq.gz"),
        html=report(
            "{species}/results/reads/reads_trimmed/fastp_report/{sample}_trimmed.pe.html",
            #caption="../report/fastp.rst",
            category="quality control",
            subcategory="fastp",
        ),
        json="{species}/results/reads/reads_trimmed/fastp_report/{sample}_trimmed.pe.json",
    message: "Trimming adapters from paired-end reads and merging for {input.sample}"
    log:
        "{species}/logs/reads/reads_trimmed/{sample}_trimmed.pe.log",
    params:
        adapters=f"--adapter_sequence {config['pipeline']['raw_reads_processing']['adapter_removal']['settings']['adapters_sequences']['r1']} --adapter_sequence_r2 {config['pipeline']['raw_reads_processing']['adapter_removal']['settings']['adapters_sequences']['r2']}",
        extra=f"--length_required {config['pipeline']['raw_reads_processing']['adapter_removal']['settings'].get('min_length','15')} --trim_poly_x 5 --qualified_quality_phred {config['pipeline']['raw_reads_processing']['adapter_removal']['settings'].get('min_quality','5')} --unqualified_percent_limit 40 --n_base_limit 5 --merge",
    threads: workflow.cores
    wrapper:
        "v7.5.0/bio/fastp"