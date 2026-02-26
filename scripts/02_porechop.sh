#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# ============================================================
# Step 2: Adapter trimming with Porechop
# Submit as: sbatch --array=1-66 --dependency=afterok:<JOB1_ID> 02_porechop.sh
# ============================================================
INPUT_DIR="/path/to/project/processing/samtools"
OUTPUT_DIR="/path/to/project/processing/porechop"
FILELIST="/path/to/project/filelist.txt"
mkdir -p "$OUTPUT_DIR"
# Get the filename for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(basename "$BAM_FILE" .bam)
INPUT_FASTQ="${INPUT_DIR}/${BASENAME}.fastq"
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "ERROR: Input file $INPUT_FASTQ does not exist. Skipping."
    exit 1
fi
echo "Processing: ${BASENAME}"
echo "Start time: $(date)"
porechop -i "$INPUT_FASTQ" \
         -o "${OUTPUT_DIR}/trimmed_${BASENAME}.fastq" \
         -t 20
echo "Finished: $(date)"
echo "Output size: $(du -h ${OUTPUT_DIR}/trimmed_${BASENAME}.fastq)"
