
# ==============================================================================
# Rule for collecting library insert size metrics using GATK/Picard
#
# Analyzes paired-end alignments to compute insert size distribution metrics
# and generates a PDF histogram plot of the results.
# ==============================================================================

rule InsertSize:
    input:
        bam=wrkdir / "alignments" / "{sample}_dedup.recall.cram",
        bai=wrkdir / "alignments" / "{sample}_dedup.recall.cram.crai",
        genome=genome,
    output:
        size_metric=wrkdir / "metrics" / "{sample}_insert_size_metrics.txt",
        pdf=wrkdir / "metrics" / "{sample}_insert_size.pdf",
    conda:
        "../envs/gatk.yaml"
    threads: 1
    resources:
        mem_mb=20000,
        runtime=24 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    log:
        logdir / "picard/{sample}_insert_size.log",
    message:
        "Collecting insert size metrics"
    shell:
        "gatk CollectInsertSizeMetrics -I {input.bam} -O {output.size_metric} -H {output.pdf} -R {input.genome} &> {log}"
