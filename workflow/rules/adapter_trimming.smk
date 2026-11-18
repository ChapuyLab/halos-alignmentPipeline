# ==============================================================================
# Rules for FASTQ preparation (BAM / SRA), adapter trimming (Cutadapt), and FastQC
# ==============================================================================


rule create_adapter_fastq:
    params:
        adapt_1=adapter_seq_r1,
        adapt_2=adapter_seq_r2,
    output:
        adapt_1=temp(scratch_path / "{sample}" / "cutadapt" / "adapt_1.fastq"),
        adapt_2=temp(scratch_path / "{sample}" / "cutadapt" / "adapt_2.fastq"),
    threads: 1
    resources:
        mem_mb=1000,
        runtime=20,
        nodes=1,
        tmpdir=scratch_dir,
    message:
        "Creating adapter fastq files"
    run:
        with open(output.adapt_1, "w") as handle:
            count = 1
            for i in params.adapt_1:
                handle.write(">adapter_" + str(count) + "\n")
                handle.write(i + "\n")
                count += 1

        with open(output.adapt_2, "w") as handle:
            count = 1
            for i in params.adapt_2:
                handle.write(">adapter_" + str(count) + "\n")
                handle.write(i + "\n")
                count += 1


rule prepare_fastq:
    input:
        bam=bamfile if input_mode == "bam" else [],
    params:
        script=lambda wc: str(Path(workflow.basedir) / "scripts" / "prepare_fastq.py"),
        mode=input_mode,
        bam=bamfile or "",
        samtools_params=extra_samtools_convert_params,
        sra_id=sra_id or "",
        sra_workdir=lambda wc: str(wrkdir / "fastq" / "sra_tmp" / wc.sample),
        ngc=ngc_file or "",
    output:
        fastq_r1=temp(scratch_path / "fastq" / "raw" / "{run_id}" / "{sample}_R1_{lane}.fastq.gz"),
        fastq_r2=temp(scratch_path / "fastq" / "raw" / "{run_id}" / "{sample}_R2_{lane}.fastq.gz"),
    log:
        logdir / "prepare_fastq/{run_id}_{sample}_{lane}.log",
    threads: 8
    resources:
        mem_mb=100000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/sra_tools.yaml"
    message:
        "Preparing paired-end FASTQ files from input source"
    shell:
        """
        python {params.script} \
            --mode {params.mode} \
            --r1 {output.fastq_r1} \
            --r2 {output.fastq_r2} \
            --threads {threads} \
            --bam "{params.bam}" \
            --samtools-params "{params.samtools_params}" \
            --sra "{params.sra_id}" \
            --sra-workdir "{params.sra_workdir}" \
            --ngc "{params.ngc}" &> {log}
        """


rule cutadapt:
    input:
        adapt_1=scratch_path / "{sample}" / "cutadapt" / "adapt_1.fastq",
        adapt_2=scratch_path / "{sample}" / "cutadapt" / "adapt_2.fastq",
        fastq_r1=lambda wc: getFastq(wc, metadata, scratch_path, False)[0],
        fastq_r2=lambda wc: getFastq(wc, metadata, scratch_path, False)[1],
    params:
        extra_cutadapt_params=extra_cutadapt_params,
    output:
        fastq_r1=temp(
            scratch_path
            / "fastq"
            / "{run_id}"
            / "cutadapt"
            / "{sample}_R1_{lane}_trim.fastq.gz"
        ),
        fastq_r2=temp(
            scratch_path
            / "fastq"
            / "{run_id}"
            / "cutadapt"
            / "{sample}_R2_{lane}_trim.fastq.gz"
        ),
    log:
        logdir / "cutadapt/{run_id}_{sample}_{lane}.log",
    threads: 8
    resources:
        mem_mb=8000,
        runtime=72 * 60,
        nodes=1,
        tmpdir=scratch_dir,
    conda:
        "../envs/cutadapt.yaml"
    message:
        "Trimming adapters using cutadapt"
    shell:
        "cutadapt {params.extra_cutadapt_params} -j {threads} -a file:{input.adapt_1} -A file:{input.adapt_2} -o {output.fastq_r1} -p {output.fastq_r2} {input.fastq_r1} {input.fastq_r2} &> {log}"


rule fastqc:
    input:
        fastq=(
            scratch_path
            / "fastq"
            / "{run_id}"
            / "cutadapt"
            / "{sample}_{read}_{lane}_trim.fastq.gz"
        ),
    output:
        html_r1=(
            wrkdir
            / "fastq"
            / "{run_id}"
            / "cutadapt"
            / "{sample}_{read}_{lane}_trim_fastqc.html"
        ),
        zip_r1=(
            wrkdir
            / "fastq"
            / "{run_id}"
            / "cutadapt"
            / "{sample}_{read}_{lane}_trim_fastqc.zip"
        ),
    log:
        logdir / "fastqc/{run_id}_{sample}_{read}_{lane}.log",
    threads: 1
    resources:
        mem_mb=1000,
        runtime=lambda wildcards, attempt: 60 * attempt,
        nodes=1,
    conda:
        "../envs/fastqc.yaml"
    message:
        "Running Fastqc on adapter trimmed sequence"
    shell:
        "fastqc {input.fastq} &> {log}"
