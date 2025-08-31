def get_reads(wc, readpair="R1"):
    """
    Returns the correct file for R1/R2, with optional suffix (_001, etc).
    """
    base = f"{wc.species}/raw/reads/{wc.individual}_{wc.rest}_{readpair}"

    # look for files that match base + (maybe something) + .fastq.gz
    candidates = [f for f in os.listdir(os.path.dirname(base))
                  if re.match(os.path.basename(base) + r"(\S*)\.fastq\.gz", f)]
        
    if not candidates:
        return None
    # pick the first candidate (or implement custom logic if multiple)
    return os.path.join(os.path.dirname(base), sorted(candidates)[0])

rule adapter_removal:
    input:
        r1 = lambda wc: str(get_reads(wc, "R1")),
        r2 = lambda wc: str(get_reads(wc, "R2"))
    output:
        trimmed="{species}/processed/trimmed/{individual}_{rest}_trimmed.fastq.gz",
        report_json="{species}/processed/trimmed/{individual}_{rest}_report.json",
        report_html="{species}/processed/trimmed/{individual}_{rest}_report.html"
    threads: workflow.cores
    params: 
        adapter_sequence_r1 = config["pipeline"]["raw_reads_processing"]["adapter_removal"]["settings"]["adapters_sequences"]["r1"],
        adapter_sequence_r2 = config["pipeline"]["raw_reads_processing"]["adapter_removal"]["settings"]["adapters_sequences"]["r2"]
    conda:
        "../../../envs/fastp.yaml"
    shell:
        """
        if [ "{input.r2}" != "None" ]; then
            fastp \
                --adapter_sequence {params.adapter_sequence_r1} \
                --adapter_sequence_r2 {params.adapter_sequence_r2} \
                --in1 {input.r1} --in2 {input.r2} \
                --merged_out {output.trimmed} \
                --json {output.report_json} \
                --html {output.report_html} \
                --merge \
                --thread {threads} \
                --length_required 15 \
                --trim_poly_x 5 \
                --qualified_quality_phred 5 \
                --unqualified_percent_limit 40 \
                --n_base_limit 5
        else
            fastp \
                --adapter_sequence {params.adapter_sequence_r1} \
                -i {input.r1} \
                -o {output.trimmed} \
                --json {output.report_json} \
                --html {output.report_html} \
                --thread {threads} \
                --length_required 15 \
                --trim_poly_x 5 \
                --qualified_quality_phred 5 \
                --unqualified_percent_limit 40 \
                --n_base_limit 5
        fi
        """
