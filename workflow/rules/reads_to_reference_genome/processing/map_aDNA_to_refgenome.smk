import scripts.setup_snakemake as setup

# 1 Map reads to reference (SAM output)
rule bwa_map_aDNA:
    input:
        reads="{species}/processed/merged/{individual}.fastq.gz",
        ref="{species}/raw/ref_genome/{ref_genome}.fa",
        index="{species}/raw/ref_genome/{ref_genome}.fa.bwt"
    output:
        sam=temp("{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}.sam")
    threads: workflow.cores 
    conda:
        os.path.join(setup.project_root, "envs/bwa.yaml")
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.reads} > {output.sam}
        """

# 2 Convert SAM to BAM
rule sam_to_bam:
    input:
        sam="{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}.sam"
    output:
        bam=temp("{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}.bam")
    threads: workflow.cores 
    conda:
        os.path.join(setup.project_root, "envs/samtools.yaml")
    shell:
        """
        samtools view -@ {threads} -bS {input.sam} -o {output.bam}
        """

# 3 Sort BAM
rule sort_bam:
    input:
        bam="{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}.bam"
    output:
        sorted_bam="{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.bam"
    threads: workflow.cores 
    conda:
        os.path.join(setup.project_root, "envs/samtools.yaml")
    shell:
        """
        samtools sort -@ {threads} {input.bam} -o {output.sorted_bam}
        """

# 4 Index BAM
rule index_bam:
    input:
        sorted_bam="{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.bam"
    output:
        bam_index="{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.bam.bai"
    threads: workflow.cores 
    conda:
        os.path.join(setup.project_root, "envs/samtools.yaml")
    shell:
        """
        samtools index -@ {threads} {input.sorted_bam} {output.bam_index}
        """
