#!/bin/bash
# -----------------------------------------------------------------------------
# Step 3b: Build read-length distributions for raw and filtered reads.
# Input: processing/samtools/<sample>.fastq and processing/nanofilt/filtered_<sample>.fastq
# Output: processing/read_lengths_raw_agg.tsv and processing/read_lengths_filtered_agg.tsv
# Run: sbatch --array=1-N --dependency=afterok:<NANOFILT_JOB> scripts/03b_read_lengths.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

RAW_DIR="${WORK_DIR}/processing/samtools"
FILT_DIR="${WORK_DIR}/processing/nanofilt"
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(derive_basename "$BAM_FILE")
RAW_FASTQ="${RAW_DIR}/${BASENAME}.fastq"
FILT_FASTQ="${FILT_DIR}/filtered_${BASENAME}.fastq"
echo "Processing: ${BASENAME}"
echo "Start time: $(date)"
# Raw read lengths (aggregated: sample \t length \t state \t count)
if [ -f "$RAW_FASTQ" ]; then
    RAW_OUT="${WORK_DIR}/processing/read_lengths_raw_agg.tsv"
    TMP_RAW=$(mktemp)
    awk -v sample="$BASENAME" 'NR%4==2 {lens[length($0)]++} END {for (l in lens) print sample "\t" l "\t" "raw" "\t" lens[l]}' \
        "$RAW_FASTQ" > "$TMP_RAW"
    (
        flock -x 200
        cat "$TMP_RAW" >> "$RAW_OUT"
    ) 200>"${RAW_OUT}.lock"
    rm -f "$TMP_RAW"
    echo "Raw lengths written"
else
    echo "WARNING: Raw FASTQ not found: $RAW_FASTQ"
fi
# Filtered read lengths (aggregated)
if [ -f "$FILT_FASTQ" ]; then
    FILT_OUT="${WORK_DIR}/processing/read_lengths_filtered_agg.tsv"
    TMP_FILT=$(mktemp)
    awk -v sample="$BASENAME" 'NR%4==2 {lens[length($0)]++} END {for (l in lens) print sample "\t" l "\t" "filtered" "\t" lens[l]}' \
        "$FILT_FASTQ" > "$TMP_FILT"
    (
        flock -x 200
        cat "$TMP_FILT" >> "$FILT_OUT"
    ) 200>"${FILT_OUT}.lock"
    rm -f "$TMP_FILT"
    echo "Filtered lengths written"
else
    echo "WARNING: Filtered FASTQ not found: $FILT_FASTQ"
fi
echo "Finished: $(date)"
