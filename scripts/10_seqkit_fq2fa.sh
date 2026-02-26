#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# -----------------------------------------------------------------------------
# Step 10: Convert filtered reads from FASTQ to FASTA format.
# Input: processing/nanofilt/filtered_<sample>.fastq
# Output: processing/fasta/<sample>.fasta
# Run: sbatch --array=1-N scripts/10_seqkit_fq2fa.sh
# -----------------------------------------------------------------------------
INPUT_DIR="/path/to/project/processing/nanofilt"
OUTPUT_DIR="/path/to/project/processing/fasta"
FILELIST="/path/to/project/filelist.txt"
mkdir -p "$OUTPUT_DIR"
# Get the filename for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(basename "$BAM_FILE" .bam)
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
