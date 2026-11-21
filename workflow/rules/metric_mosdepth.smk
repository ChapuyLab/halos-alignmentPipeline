# ==============================================================================
# Rules for depth & coverage metric calculation using Mosdepth
#
# Runs mosdepth either in WGS mode (no target regions) or Panel/WES mode
# (using the target BED file to restrict calculation to targeted intervals).
# ==============================================================================

rule mosdepth:
    input:
        bam=wrkdir / "alignments" / "{sample}_dedup.recall.cram",
        bai=wrkdir / "alignments" / "{sample}_dedup.recall.cram.crai",
        genome=genome,
        target_regions=mosdepth_target_regions,
    params:
        prefix=str(wrkdir / "metrics" / "{sample}"),
        by_flag=mosdepth_by_flag,
    output:
        out_1=wrkdir / "metrics" / "{sample}.mosdepth.global.dist.txt",
        out_2=wrkdir / "metrics" / "{sample}.mosdepth.summary.txt",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=24 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/mosdepth.yaml"
    log:
        logdir / "mosdepth/{sample}.log",
    message:
        "Running mosdepth for WES/panel data" if mosdepth_by_flag else "Running mosdepth for WGS data"
    shell:
        "mosdepth -f {input.genome} {params.by_flag} -n {params.prefix} {input.bam} &> {log}"
