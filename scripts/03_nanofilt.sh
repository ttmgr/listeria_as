#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# ============================================================
# Step 3: Length filtering with NanoFilt (remove reads < 100 bp)
# Submit as: sbatch --array=1-66 --dependency=afterok:<JOB2_ID> 03_nanofilt.sh
# ============================================================
INPUT_DIR="/path/to/project/processing/porechop"
OUTPUT_DIR="/path/to/project/processing/nanofilt"
FILELIST="/path/to/project/filelist.txt"
mkdir -p "$OUTPUT_DIR"
# Get the filename for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(basename "$BAM_FILE" .bam)
INPUT_FASTQ="${INPUT_DIR}/trimmed_${BASENAME}.fastq"
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "ERROR: Input file $INPUT_FASTQ does not exist. Skipping."
    exit 1
fi
echo "Processing: ${BASENAME}"
echo "Start time: $(date)"
cat "$INPUT_FASTQ" | NanoFilt -l 100 > "${OUTPUT_DIR}/filtered_${BASENAME}.fastq"
# Quick stats
INPUT_READS=$(echo $(cat "$INPUT_FASTQ" | wc -l) / 4 | bc)
OUTPUT_READS=$(echo $(cat "${OUTPUT_DIR}/filtered_${BASENAME}.fastq" | wc -l) / 4 | bc)
echo "Reads in:  $INPUT_READS"
echo "Reads out: $OUTPUT_READS"
echo "Removed:   $((INPUT_READS - OUTPUT_READS)) reads below 100 bp"
echo "Finished: $(date)"
