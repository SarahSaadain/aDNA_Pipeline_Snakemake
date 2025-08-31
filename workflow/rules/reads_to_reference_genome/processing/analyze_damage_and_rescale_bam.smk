import scripts.setup_snakemake as setup
rule analyze_damage_and_rescale_bam:
    """
    Estimates post-mortem DNA damage patterns and rescales base qualities.
    """
    input:
        bam = "{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.bam",
        bam_index = "{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.bam.bai",
        ref = "{species}/raw/ref_genome/{ref_genome}.fa"
    output:
        log = "{species}/results/{ref_genome}/damage/{individual}/Runtime_log.txt",
        GtoA3p = "{species}/results/{ref_genome}/damage/{individual}/3pGtoA_freq.txt",
        CtoT5p = "{species}/results/{ref_genome}/damage/{individual}/5pCtoT_freq.txt",
        dnacomp = "{species}/results/{ref_genome}/damage/{individual}/dnacomp.txt",
        frag_misincorp = "{species}/results/{ref_genome}/damage/{individual}/Fragmisincorporation_plot.pdf",
        #len = "{species}/results/{ref_genome}/damage/{individual}/Length_plot.pdf",
        lg_dist = "{species}/results/{ref_genome}/damage/{individual}/lgdistribution.txt",
        misincorp = "{species}/results/{ref_genome}/damage/{individual}/misincorporation.txt",
        rescaled_bam = temp("{species}/results/{ref_genome}/damage/{individual}/{individual}_{ref_genome}.bam")
    params:
        extra="--merge-reference-sequences --rescale",  # optional parameters for mapdamage2 (except -i, -r, -d, --rescale)
    message:
        "Analyze damage and rescale {input.bam}",
    benchmark:
        "benchmark/{species}_{individual}_{ref_genome}_mapdamage2.benchmark.json",
    log:
        "logs/{species}/{individual}_{ref_genome}_mapdamage2.log",
    wrapper:
        "v7.2.0/bio/mapdamage2"

# 3 Sort BAM
rule sort_rescaled_bam:
    input:
        rescaled_bam="{species}/results/{ref_genome}/damage/{individual}/{individual}_{ref_genome}.bam"
    output:
        sorted_bam="{species}/results/{ref_genome}/damage/{individual}/{individual}_{ref_genome}_sorted.bam"
    threads: workflow.cores 
    conda:
        os.path.join(setup.project_root, "envs/samtools.yaml")
    shell:
        """
        samtools sort -@ {threads} {input.rescaled_bam} -o {output.sorted_bam}
        """

# 4 Index BAM
rule index_rescaled_bam:
    input:
        sorted_bam="{species}/results/{ref_genome}/damage/{individual}/{individual}_{ref_genome}_sorted.bam"
    output:
        bam_index="{species}/results/{ref_genome}/damage/{individual}/{individual}_{ref_genome}_sorted.bam.bai"
    threads: workflow.cores 
    conda:
        os.path.join(setup.project_root, "envs/samtools.yaml")
    shell:
        """
        samtools index -@ {threads} {input.sorted_bam} {output.bam_index}
        """

rule move_rescaled_bam:
    input:
        sorted_bam="{species}/results/{ref_genome}/damage/{individual}/{individual}_{ref_genome}_sorted.bam",
        bam_index="{species}/results/{ref_genome}/damage/{individual}/{individual}_{ref_genome}_sorted.bam.bai"
    output:
        sorted_bam="{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.rescaled.bam",
        bam_index="{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.rescaled.bam.bai"
    shell:
        """
        mv {input.sorted_bam} {output.sorted_bam} 
        mv {input.bam_index} {output.bam_index} 
        """

