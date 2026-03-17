#!/bin/bash
# -----------------------------------------------------------------------------
# Step 10: Convert filtered reads from FASTQ to FASTA format.
# Input: processing/nanofilt/filtered_<sample>.fastq
# Output: processing/fasta/<sample>.fasta
# Run: sbatch --array=1-N scripts/10_seqkit_fq2fa.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

INPUT_DIR="${WORK_DIR}/processing/nanofilt"
OUTPUT_DIR="${WORK_DIR}/processing/fasta"
mkdir -p "$OUTPUT_DIR"
# Get the filename for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(derive_basename "$BAM_FILE")
INPUT_FASTQ="${INPUT_DIR}/filtered_${BASENAME}.fastq"
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "ERROR: Input file $INPUT_FASTQ does not exist. Skipping."
    exit 1
fi
echo "Processing: ${BASENAME}"
echo "Start time: $(date)"
seqkit fq2fa "$INPUT_FASTQ" -o "${OUTPUT_DIR}/${BASENAME}.fasta"
echo "Finished: $(date)"
echo "Output size: $(du -h ${OUTPUT_DIR}/${BASENAME}.fasta)"
