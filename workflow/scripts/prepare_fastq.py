#!/usr/bin/env python3
"""Prepare paired-end FASTQ files from either a BAM file or an NCBI SRA accession."""

import argparse
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path


def run_bam_to_fastq(bam: str, r1: Path, r2: Path, threads: int, samtools_params: str) -> None:
    """Convert an aligned/unaligned BAM file to paired-end gzipped FASTQs."""
    if not bam:
        raise ValueError("--bam is required when --mode is 'bam'")

    extra_args = shlex.split(samtools_params) if samtools_params else []
    with tempfile.TemporaryDirectory(prefix="bam_collate_", dir=str(r1.parent)) as collate_dir:
        reheader_cmd = ["samtools", "reheader", "-c", "sed 's/MT5:/M5:/g'", bam]
        collate_cmd = ["samtools", "collate", "-T", str(Path(collate_dir) / "collate"), "-u", "-O", "-@", str(threads), "-"]
        fastq_cmd = [
            "samtools", "fastq", *extra_args,
            "-n", "-F", "0x900",
            "-1", str(r1), "-2", str(r2),
            "-0", "/dev/null", "-s", "/dev/null", "-",
        ]

        p1 = subprocess.Popen(reheader_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(collate_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        p3 = subprocess.Popen(fastq_cmd, stdin=p2.stdout)
        p2.stdout.close()

        rc3 = p3.wait()
        rc2 = p2.wait()
        rc1 = p1.wait()
        if rc1 != 0 or rc2 != 0 or rc3 != 0:
            raise RuntimeError(f"BAM-to-FASTQ pipeline failed with exit codes: reheader={rc1}, collate={rc2}, fastq={rc3}")


def run_sra_to_fastq(sra_id: str, r1: Path, r2: Path, threads: int, sra_workdir: Path, ngc: str) -> None:
    """Download an SRA accession and convert it to paired-end gzipped FASTQs inside sra_workdir."""
    if not sra_id:
        raise ValueError("--sra is required when --mode is 'sra'")

    ngc_args = ["--ngc", ngc] if ngc else []
    sra_workdir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f"sra_{sra_id}_", dir=str(sra_workdir)) as sra_tmp:
        sra_tmp_path = Path(sra_tmp)
        prefetch_cmd = ["prefetch", *ngc_args, sra_id, "--max-size", "u", "-O", str(sra_tmp_path)]
        subprocess.run(prefetch_cmd, check=True)

        fasterq_cmd = [
            "fasterq-dump", *ngc_args,
            "--split-files",
            "--threads", str(threads),
            "-t", str(sra_tmp_path),
            "-O", str(sra_tmp_path),
            str(sra_tmp_path / sra_id),
        ]
        subprocess.run(fasterq_cmd, check=True)

        fq1 = sra_tmp_path / f"{sra_id}_1.fastq"
        fq2 = sra_tmp_path / f"{sra_id}_2.fastq"
        if not fq1.is_file() or not fq2.is_file():
            raise FileNotFoundError(
                f"Expected paired-end FASTQs ({fq1.name} and {fq2.name}) were not generated for {sra_id}."
            )

        for src, dst in [(fq1, r1), (fq2, r2)]:
            with open(dst, "wb") as out_f:
                subprocess.run(["pigz", "-p", str(threads), "-c", str(src)], stdout=out_f, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare paired-end FASTQ files from BAM or SRA.")
    parser.add_argument("--mode", required=True, choices=["bam", "sra"], help="Input mode: 'bam' or 'sra'")
    parser.add_argument("--r1", required=True, help="Output R1 .fastq.gz path")
    parser.add_argument("--r2", required=True, help="Output R2 .fastq.gz path")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    parser.add_argument("--bam", default="", help="Input BAM path (for mode='bam')")
    parser.add_argument("--samtools-params", default="", help="Extra flags for samtools fastq")
    parser.add_argument("--sra", default="", help="SRA accession (for mode='sra')")
    parser.add_argument("--sra-workdir", default="", help="Working directory in work_dir for SRA download/dump (for mode='sra')")
    parser.add_argument("--ngc", default="", help="Optional .ngc key file (for mode='sra')")

    args = parser.parse_args()
    r1 = Path(args.r1)
    r2 = Path(args.r2)
    sra_workdir = Path(args.sra_workdir) if args.sra_workdir else (r1.parent / "sra_tmp")

    r1.parent.mkdir(parents=True, exist_ok=True)
    r2.parent.mkdir(parents=True, exist_ok=True)

    if args.mode == "bam":
        run_bam_to_fastq(args.bam, r1, r2, args.threads, args.samtools_params)
    else:
        run_sra_to_fastq(args.sra, r1, r2, args.threads, sra_workdir, args.ngc)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(1)
