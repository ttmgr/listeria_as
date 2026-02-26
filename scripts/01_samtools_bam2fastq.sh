#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# ============================================================
# Step 1: BAM → FASTQ conversion using samtools
# Submit as: sbatch --array=1-66 01_samtools_bam2fastq.sh
# ============================================================
INPUT_DIR="/path/to/BAM_files"
OUTPUT_DIR="/path/to/project/processing/samtools"
FILELIST="/path/to/project/filelist.txt"
mkdir -p "$OUTPUT_DIR"
# Get the BAM file for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
if [ -z "$BAM_FILE" ]; then
    echo "ERROR: No file found for task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi
BASENAME=$(basename "$BAM_FILE" .bam)
echo "Processing: ${BAM_FILE} → ${BASENAME}.fastq"
echo "Start time: $(date)"
samtools fastq -@ 4 "${INPUT_DIR}/${BAM_FILE}" > "${OUTPUT_DIR}/${BASENAME}.fastq"
echo "Finished: $(date)"
echo "Output size: $(du -h ${OUTPUT_DIR}/${BASENAME}.fastq)"
