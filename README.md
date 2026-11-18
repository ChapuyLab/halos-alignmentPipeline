# Alignment Pipeline 2.0 (FASTQ / BAM / SRA)

Author: Shashwat Sahay (shashwat.sahay@charite.de)

Unified paired-end alignment, realignment, and QC pipeline for **WGS**, **WES**, and **Panel** sequencing data. Supports three input modes from a single Snakemake workflow:

1. **FASTQ metadata sheet (`metadata`)** — multi-run / multi-lane paired-end FASTQs
2. **BAM file (`bamfile`)** — extracts paired-end reads and realigns from scratch
3. **SRA accession (`sra` / `sra_id`)** — downloads and converts SRA runs to paired-end FASTQs (with optional `.ngc` key file for controlled-access data)

> **Requirement:** Sequencing data must be paired-end.

---

## Prerequisites

- A working **Snakemake** installation in a conda/mamba environment (recommended: `miniforge3`).
- Per-rule conda environments are defined in `workflow/envs/` and managed automatically by Snakemake when run with `--use-conda`.
  - `workflow/envs/sra_tools.yaml` provides `sra-tools`, `samtools`, and `pigz` for BAM/SRA FASTQ preparation.

---

## Input Modes

Provide **exactly one** of `metadata`, `bamfile`, or `sra` (or `sra_id`) in `config/config.yaml`:

- **FASTQ Metadata (`metadata`)**: Path to a CSV/TSV metadata sheet with columns `PATIENT_ID`, `SAMPLE_TYPE`, `SAMPLE_NAME`, `RUN_ID`, `LANE_NO`, `READ` (`R1`/`R2`), and `FASTQ_FILE`. Supports multiple runs and lanes, aligned per lane and merged before duplicate marking.
- **BAM (`bamfile`)**: Path to an aligned/unaligned BAM file. `workflow/scripts/prepare_fastq.sh` fixes `MT5:` headers if present, collates pairs, and extracts paired-end FASTQs via `samtools fastq`.
- **SRA (`sra` or `sra_id`)**: NCBI SRA run accession (e.g., `SRR12345678`). `workflow/scripts/prepare_fastq.sh` downloads the run with `prefetch`, converts it to paired-end FASTQs with `fasterq-dump --split-files`, and compresses with `pigz`.
  - **Controlled-access SRA (`ngc` or `ngc_file`, optional)**: Path to a dbGaP/SRA `.ngc` key file passed via `--ngc` to `prefetch` and `fasterq-dump`.

---

## Configuration (`config/config.yaml`)

- **Sample & Sequencing**: `SeqType` (`"WGS"`, `"WES"`, `"Panel"`), `library_prep_kit`, `pid`, `sample`.
- **Input**: One of `metadata`, `bamfile`, or `sra` (plus optional `ngc` for SRA).
- **Directories**: `work_dir`, `log_dir` (defaults to `{work_dir}/logs`), `scratch_dir` (optional; intermediate FASTQs/BAMs are stored here, defaulting to `work_dir`).
- **Reference & Targets**: `genome`, `dbsnp`, and (for `WES`/`Panel`) `target_regions`, `bait_regions` (optional), `chrom_sizes`, `dict_genome`.
- **Trimming & Options**: `trim_adapters` (`True`/`False`), `adapter_forward_read`, `adapter_reverse_read`, `extra_cutadapt_params`, `extra_samtools_convert_params`, `max_coverage`.
  > **Note (`extra_cutadapt_params`):** Modern exome library preparations such as **Agilent SureSelect V8** include UMI tags at the start of the reads. If your data contains these UMI tags, set `extra_cutadapt_params: "-u 5"` (with `trim_adapters: True`) to trim them off before alignment.

---

## Running the Pipeline

From the `alignment_pipeline_2.0` directory:

```bash
snakemake --use-conda --cores <N> --configfile config/config.yaml
```

Or using the Slurm profile:

```bash
snakemake --profile profile --configfile config/config.yaml
```

---

## Outputs

- **Final Alignment (`{work_dir}/alignments/`)**:
  - `{pid}_{sample}_dedup.recall.cram` and `.crai`
- **QC Metrics (`{work_dir}/metrics/`)**:
  - `{pid}_{sample}.mosdepth.global.dist.txt`, `{pid}_{sample}.mosdepth.summary.txt`
  - `{pid}_{sample}_dedup.recall.flagstat`
  - `{pid}_{sample}_insert_size_metrics.txt`, `{pid}_{sample}_insert_size.pdf`
  - `{pid}_{sample}_recal_data.table`, `{pid}_{sample}_covariates.pdf`
  - `{pid}_{sample}.hs_metrics.txt` (`WES` / `Panel`) or `{pid}_{sample}_coverage.png` (`WGS`)
  - FastQC reports in `{work_dir}/fastq/{run_id}/cutadapt/` when `trim_adapters: True`
