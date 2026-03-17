#!/bin/bash
# -----------------------------------------------------------------------------
# Step 1: Convert Nanopore BAM files to FASTQ.
# Input: BAM file chosen by SLURM array index from filelist.txt
# Output: processing/samtools/<sample>.fastq
# Run: sbatch --array=1-N scripts/01_samtools_bam2fastq.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

OUTPUT_DIR="${WORK_DIR}/processing/samtools"
mkdir -p "$OUTPUT_DIR"
# Get the BAM file for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
if [ -z "$BAM_FILE" ]; then
    echo "ERROR: No file found for task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi
BASENAME=$(derive_basename "$BAM_FILE")
echo "Processing: ${BAM_FILE} → ${BASENAME}.fastq"
echo "Start time: $(date)"
samtools fastq -@ 4 "${BAM_BASE_DIR}/${BAM_FILE}" > "${OUTPUT_DIR}/${BASENAME}.fastq"
echo "Finished: $(date)"
echo "Output size: $(du -h ${OUTPUT_DIR}/${BASENAME}.fastq)"
