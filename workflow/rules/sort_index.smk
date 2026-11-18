# ==============================================================================
# Rule for sorting and indexing BAM files using Samtools
#
# Coordinate-sorts BAM alignments and creates a corresponding BAI index file.
# ==============================================================================

rule sort_index:
    input:
        bam=scratch_path / "alignments" / "{sample}.bam",
    output:
        bam=temp(scratch_path / "alignments" / "{sample}.sorted.bam"),
        bai=temp(scratch_path / "alignments" / "{sample}.sorted.bam.bai"),
    conda:
        "../envs/samtools.yaml"
    threads: 8
    params:
        mem_thread=8000,
    resources:
        mem_mb=lambda wildcards, threads: (threads + 1) * 8000,
        runtime=24 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        logdir / "samtools/{sample}_sort.log",
    message:
        "Sorting and indexing recalibrated bam file"
    shell:
        " ( "
        " samtools sort --threads {threads} -m{params.mem_thread}m -o {output.bam}##idx##{output.bai} {input.bam} -T {resources.tmpdir} --write-index"
        " ) &> {log} "
