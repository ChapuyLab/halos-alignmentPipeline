# ==============================================================================
# Rules for alignment and lane-merging using BWA and Samtools
# ==============================================================================


rule bwa_map:
    """
    First pass alignemnt
    Aligning reads to the genome using BWA
    """
    input:
        genome=genome,
        fastq_r1=lambda wc: getFastq(wc, metadata, scratch_path, config["trim_adapters"]),
    params:
        RG=lambda wc: "'@RG\\tID:"
        + wc.run_id
        + "_"
        + wc.lane
        + "\\tSM:"
        + config["pid"]
        + "_"
        + config["sample"]
        + "\\tLB:"
        + config["pid"]
        + "_"
        + config["sample"]
        + "\\tPL:ILLUMINA'",
    output:
        temp(scratch_path / "alignments" / "{run_id}" / "{sample}_aln_{lane}.bam"),
    threads: 8
    resources:
        mem_mb=20000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/samtools.yaml"
    message:
        "First pass alignemnt. Aligning reads to the genome using BWA."
    log:
        logdir / "bwa" / "first_pass_align_{run_id}_{sample}_{lane}.log",
    shell:
        "bwa mem -Y -K 150000000 -t {threads} -R {params.RG} {input.genome} {input.fastq_r1} "
        "| samtools view -b - > {output}"


rule merge:
    """
    Merging bam files from different lanes/runs
    """
    input:
        expand(
            scratch_path / "alignments" / "{run_id}" / "{sample}_aln_{lane}.bam",
            filtered_product,
            run_id=RUN_ID,
            sample=config["pid"] + "_" + config["sample"],
            lane=LANE,
        ),
    output:
        bam=temp(scratch_path / "alignments" / "{sample}.merged.bam"),
    threads: 24
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/samtools.yaml"
    message:
        "Merging bam files from different lanes/runs."
    log:
        logdir / "samtools/{sample}_merge.log",
    shell:
        " if [ $(echo {input} | wc -w) -eq 1 ]; then "
        " mv {input[0]} {output.bam}; "
        " else "
        "  samtools merge --threads {threads} -f {output.bam} {input}; "
        "fi"
