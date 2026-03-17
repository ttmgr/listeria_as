#!/bin/bash
# -----------------------------------------------------------------------------
# Step 5: Classify each filtered read taxonomically with Kraken2.
# Input: processing/nanofilt/filtered_<sample>.fastq
# Output: processing/kraken2/report_<sample>.txt and classified_<sample>.txt
# Run: sbatch --array=1-N scripts/05_kraken2.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

INPUT_DIR="${WORK_DIR}/processing/nanofilt"
OUTPUT_DIR="${WORK_DIR}/processing/kraken2"
mkdir -p "$OUTPUT_DIR"
# Get the filename for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(derive_basename "$BAM_FILE")
INPUT_FASTQ="${INPUT_DIR}/filtered_${BASENAME}.fastq"
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "ERROR: Input file $INPUT_FASTQ does not exist. Skipping."
    exit 1
fi
KRAKEN2_REPORT="${OUTPUT_DIR}/report_${BASENAME}.txt"
KRAKEN2_OUTPUT="${OUTPUT_DIR}/classified_${BASENAME}.txt"
echo "Processing: ${BASENAME}"
echo "Input: ${INPUT_FASTQ}"
echo "Start time: $(date)"
kraken2 --db "$KRAKEN2_DB" \
        --use-names \
        --threads ${SLURM_CPUS_PER_TASK:-8} \
        --report "$KRAKEN2_REPORT" \
        --output "$KRAKEN2_OUTPUT" \
        "$INPUT_FASTQ"
echo "Finished Kraken2: $(date)"
