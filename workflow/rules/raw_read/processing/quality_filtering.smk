rule quality_filter:
    input:
        trimmed="{species}/processed/trimmed/{individual}_{rest}_trimmed.fastq.gz"
    output:
        filtered="{species}/processed/quality_filtered/{individual}_{rest}_quality_filtered.fastq.gz",
        report_json="{species}/processed/quality_filtered/{individual}_{rest}_quality_filtered_report.json",
        report_html="{species}/processed/quality_filtered/{individual}_{rest}_quality_filtered_report.html",
        failed="{species}/processed/quality_filtered/{individual}_{rest}_quality_filtered_failed.fastq.gz"
    threads: workflow.cores 
    conda:
        "../../../envs/fastp.yaml"
    shell:
        """
        fastp \
            --thread {threads} \
            --disable_adapter_trimming \
            --qualified_quality_phred 15 \
            --length_required 15 \
            --unqualified_percent_limit 40 \
            --n_base_limit 5 \
            -i {input.trimmed} \
            -o {output.filtered} \
            --failed_out {output.failed} \
            --json {output.report_json} \
            --html {output.report_html}
        """
