#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# -----------------------------------------------------------------------------
# Step 4: Calculate per-sample read QC metrics (count, N50, quality, bases).
# Input: processing/nanofilt/filtered_<sample>.fastq
# Output: processing/nanostat/nanostat_<sample>.txt
# Run: sbatch --array=1-N scripts/04_nanostat.sh
# -----------------------------------------------------------------------------
INPUT_DIR="/path/to/project/processing/nanofilt"
OUTPUT_DIR="/path/to/project/processing/nanostat"
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
NanoStat --fastq "$INPUT_FASTQ" \
         --name "${BASENAME}" \
         --outdir "$OUTPUT_DIR" \
         --tsv
# Rename output to include sample name
mv "${OUTPUT_DIR}/${BASENAME}NanoStats.txt" "${OUTPUT_DIR}/nanostat_${BASENAME}.txt" 2>/dev/null || true
echo "Finished: $(date)"
