def get_reads_list(wc):
    """
    Returns a list of read files for R1/R2 if available.
    If only R1 exists, returns a single-element list [R1].
    """
    reads_dir = f"{wc.species}/raw/reads"
    
    # R1
    base_r1 = f"{wc.individual}_{wc.rest}_R1"
    candidates_r1 = [f for f in os.listdir(reads_dir)
                     if re.match(base_r1 + r"(\S*)?\.fastq\.gz", f)]
    if not candidates_r1:
        raise FileNotFoundError(f"No R1 found for {wc.individual}_{wc.rest}")
    r1 = os.path.join(reads_dir, sorted(candidates_r1)[0])
    
    # R2
    base_r2 = f"{wc.individual}_{wc.rest}_R2"
    candidates_r2 = [f for f in os.listdir(reads_dir)
                     if re.match(base_r2 + r"(\S*)?\.fastq\.gz", f)]
    
    if candidates_r2:
        r2 = os.path.join(reads_dir, sorted(candidates_r2)[0])
        return [r1, r2]  # Paired-end
    else:
        return [r1]      # Single-end

rule adapter_removal:
    input:
        reads = lambda wc: get_reads_list(wc)
    output:
        trimmed=temp("{species}/processed/trimmed/{individual}_{rest}_trimmed.fastq.gz"),
        report_json="{species}/processed/trimmed/{individual}_{rest}_report.json",
        report_html="{species}/processed/trimmed/{individual}_{rest}_report.html"
    threads: workflow.cores
    params: 
        adapter_sequence_r1 = config["pipeline"]["raw_reads_processing"]["adapter_removal"]["settings"]["adapters_sequences"]["r1"],
        adapter_sequence_r2 = config["pipeline"]["raw_reads_processing"]["adapter_removal"]["settings"]["adapters_sequences"]["r2"]
    conda:
        "../../../envs/fastp.yaml"
    run:
        reads = input.reads
        if len(reads) == 2:
            shell(f"fastp --adapter_sequence {params.adapter_sequence_r1} --adapter_sequence_r2 {params.adapter_sequence_r2} --in1 {reads[0]} --in2 {reads[1]} --merged_out {output.trimmed} --json {output.report_json} --html {output.report_html} --merge --length_required 15 --trim_poly_x 5 --qualified_quality_phred 5 --unqualified_percent_limit 40 --n_base_limit 5 --thread {threads}")
        else:
            shell(f"fastp --adapter_sequence {params.adapter_sequence_r1} -i {reads[0]} -o {output.trimmed} --json {output.report_json} --html {output.report_html} --length_required 15 --trim_poly_x 5  --qualified_quality_phred 5 --unqualified_percent_limit 40 --n_base_limit 5 --thread {threads}")