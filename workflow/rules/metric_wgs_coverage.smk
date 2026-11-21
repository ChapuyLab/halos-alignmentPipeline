# ==============================================================================
# Rules for calculating binned WGS read counts and coverage plot from CRAM
#
# 1. Counts MAPQ >= 37 primary non-duplicate read starts in 100kb bins directly
#    from CRAM using samtools selective decoding (FLAG, RNAME, POS, MAPQ only),
#    matching QDNAseq::binReadCounts() default filtering & binning behavior.
# 2. Builds a QDNAseqReadCounts RDS object and plots genome-wide coverage.
# ==============================================================================

WGS_BINSIZE = 100


rule wgsBinReadCounts:
    input:
        bam=wrkdir / "alignments" / "{sample}_dedup.recall.cram",
        bai=wrkdir / "alignments" / "{sample}_dedup.recall.cram.crai",
        genome=genome,
    params:
        bin_bp=WGS_BINSIZE * 1000,
    output:
        counts=temp(scratch_path / "metrics" / ("{sample}_" + str(WGS_BINSIZE) + "kb_read_counts.tsv")),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=12 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/samtools.yaml"
    log:
        logdir / "coveragePlot/{sample}_bin_counts.log",
    message:
        "Counting read starts in " + str(WGS_BINSIZE) + "kb bins from CRAM"
    shell:
        "("
        " set -o pipefail;"
        " samtools view -@ {threads}"
        " --input-fmt-option required_fields=0x1e"
        " -T {input.genome}"
        " -F 0x604 -q 37"
        " {input.bam}"
        " | awk -F '\\t' -v bs={params.bin_bp}"
        " 'BEGIN{{OFS=\"\\t\"}} $5 != 255 {{ b = int(($4 - 1) / bs) * bs + 1; c[$3 SUBSEP b]++ }}"
        " END{{ for (k in c) {{ split(k, a, SUBSEP); print a[1], a[2], c[k] }} }}'"
        " > {output.counts}"
        ") 2> {log}"


rule coveragePlot:
    input:
        counts=scratch_path / "metrics" / ("{sample}_" + str(WGS_BINSIZE) + "kb_read_counts.tsv"),
    params:
        binsize=WGS_BINSIZE,
    output:
        plot=wrkdir / "metrics" / "{sample}_coverage.png",
        rds=wrkdir / "metrics" / "{sample}_readCounts.rds",
    threads: 4
    resources:
        mem_mb=32000,
        runtime=24 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/coveragePlot.yaml"
    log:
        logdir / "coveragePlot/{sample}.log",
    message:
        "Building QDNAseq readCounts RDS and plotting Whole Genome coverage"
    script:
        "../scripts/getCoveragePlot_snakemake.R"
