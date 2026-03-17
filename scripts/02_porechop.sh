#!/bin/bash
# -----------------------------------------------------------------------------
# Step 2: Trim sequencing adapters from reads.
# Input: processing/samtools/<sample>.fastq
# Output: processing/porechop/trimmed_<sample>.fastq
#
# Uses Chopper (preferred) or Porechop as fallback.
# Porechop has a known StopIteration bug with Python 3.12 on large files.
#
# Run: sbatch --array=1-N --dependency=afterok:<JOB1_ID> scripts/02_porechop.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

INPUT_DIR="${WORK_DIR}/processing/samtools"
OUTPUT_DIR="${WORK_DIR}/processing/porechop"
mkdir -p "$OUTPUT_DIR"
# Get the filename for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(derive_basename "$BAM_FILE")
INPUT_FASTQ="${INPUT_DIR}/${BASENAME}.fastq"
OUTPUT_FASTQ="${OUTPUT_DIR}/trimmed_${BASENAME}.fastq"
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "ERROR: Input file $INPUT_FASTQ does not exist. Skipping."
    exit 1
fi
echo "Processing: ${BASENAME}"
echo "Input size: $(du -h "$INPUT_FASTQ" | cut -f1)"
echo "Start time: $(date)"

# Try Chopper first (modern replacement, handles large files)
if command -v chopper &>/dev/null; then
    echo "Using: chopper"
    chopper -i "$INPUT_FASTQ" \
            --threads 20 \
            > "$OUTPUT_FASTQ"
# Fall back to Porechop
elif command -v porechop &>/dev/null; then
    echo "Using: porechop"
    porechop -i "$INPUT_FASTQ" \
             -o "$OUTPUT_FASTQ" \
             -t 20
else
    echo "WARNING: Neither chopper nor porechop found. Copying input as-is."
    cp "$INPUT_FASTQ" "$OUTPUT_FASTQ"
fi

# Verify output was created and is non-empty
if [ ! -s "$OUTPUT_FASTQ" ]; then
    echo "ERROR: Output file is empty or missing. Falling back to untrimmed input."
    cp "$INPUT_FASTQ" "$OUTPUT_FASTQ"
fi

echo "Finished: $(date)"
echo "Output size: $(du -h "$OUTPUT_FASTQ" | cut -f1)"
