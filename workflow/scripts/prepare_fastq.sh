#!/usr/bin/env bash
set -euo pipefail

MODE=""
R1=""
R2=""
THREADS=8
TMPDIR="${TMPDIR:-/tmp}"
BAM=""
SAMTOOLS_PARAMS=""
SRA_ID=""
NGC=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --r1) R1="$2"; shift 2 ;;
        --r2) R2="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --tmpdir) TMPDIR="$2"; shift 2 ;;
        --bam) BAM="$2"; shift 2 ;;
        --samtools-params) SAMTOOLS_PARAMS="$2"; shift 2 ;;
        --sra) SRA_ID="$2"; shift 2 ;;
        --ngc) NGC="$2"; shift 2 ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "$(dirname "$R1")" "$(dirname "$R2")" "$TMPDIR"

if [[ "$MODE" == "bam" ]]; then
    if [[ -z "$BAM" ]]; then
        echo "[ERROR] --bam is required when --mode is 'bam'" >&2
        exit 1
    fi
    samtools reheader -c "sed 's/MT5:/M5:/g'" "$BAM" | \
        samtools collate -T "$TMPDIR" -u -O -@ "$THREADS" - | \
        samtools fastq $SAMTOOLS_PARAMS -n -F 0x900 -1 "$R1" -2 "$R2" -0 /dev/null -s /dev/null -

elif [[ "$MODE" == "sra" ]]; then
    if [[ -z "$SRA_ID" ]]; then
        echo "[ERROR] --sra is required when --mode is 'sra'" >&2
        exit 1
    fi
    SRA_TMP=$(mktemp -d "${TMPDIR}/sra_${SRA_ID}_XXXXXX")
    trap 'rm -rf "$SRA_TMP"' EXIT

    NGC_ARGS=()
    if [[ -n "$NGC" ]]; then
        NGC_ARGS=(--ngc "$NGC")
    fi

    prefetch "${NGC_ARGS[@]}" "$SRA_ID" --max-size u -O "$SRA_TMP"
    fasterq-dump "${NGC_ARGS[@]}" --split-files --threads "$THREADS" -t "$SRA_TMP" -O "$SRA_TMP" "$SRA_TMP/$SRA_ID"

    if [[ ! -f "$SRA_TMP/${SRA_ID}_1.fastq" || ! -f "$SRA_TMP/${SRA_ID}_2.fastq" ]]; then
        echo "[ERROR] Expected paired-end FASTQs (${SRA_ID}_1.fastq and _2.fastq) were not generated." >&2
        exit 1
    fi

    pigz -p "$THREADS" -c "$SRA_TMP/${SRA_ID}_1.fastq" > "$R1"
    pigz -p "$THREADS" -c "$SRA_TMP/${SRA_ID}_2.fastq" > "$R2"

else
    echo "[ERROR] Unsupported mode: '$MODE' (expected 'bam' or 'sra')" >&2
    exit 1
fi
