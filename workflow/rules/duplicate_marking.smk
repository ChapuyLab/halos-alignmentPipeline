# ==============================================================================
# Rule for duplicate marking using Sambamba
#
# Marks duplicate alignments in the lane-merged BAM file to improve BQSR
# and downstream variant calling accuracy.
# ==============================================================================

rule duplicates:
    input:
        bam=scratch_path / "alignments" / "{sample}.merged.sorted.bam",
        bai=scratch_path / "alignments" / "{sample}.merged.sorted.bam.bai",
    output:
        bam=temp(scratch_path / "alignments" / "{sample}_dedup.bam"),
    conda:
        "../envs/sambamba.yaml"
    threads: 8
    log:
        logdir / "gatk/{sample}_dedup.log",
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    message:
        "Marking duplicates"
    shell:
        " ( "
        " sambamba markdup -t {threads} --tmpdir {resources.tmpdir} {input.bam} {output.bam}"
        " ) &> {log} "
