#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# -----------------------------------------------------------------------------
# Step 8b: Assemble filtered reads with Myloasm.
# Input: processing/nanofilt/filtered_<sample>.fastq
# Output: processing/myloasm/<sample>/assembly_primary.fa
# Run: sbatch --array=1-N scripts/08b_myloasm.sh
# -----------------------------------------------------------------------------
INPUT_DIR="/path/to/project/processing/nanofilt"
OUTPUT_DIR="/path/to/project/processing/myloasm"
FILELIST="/path/to/project/filelist.txt"
# Get the filename for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(basename "$BAM_FILE" .bam)
INPUT_FASTQ="${INPUT_DIR}/filtered_${BASENAME}.fastq"
SAMPLE_OUT="${OUTPUT_DIR}/${BASENAME}"
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "ERROR: Input file $INPUT_FASTQ does not exist. Skipping."
    exit 1
fi
mkdir -p "$SAMPLE_OUT"
echo "Processing: ${BASENAME}"
echo "Start time: $(date)"
myloasm "$INPUT_FASTQ" \
    -o "$SAMPLE_OUT" \
    -t 20
# Fallback: if Myloasm aborted due to lack of reads, create a dummy file
if [ ! -f "${SAMPLE_OUT}/assembly_primary.fa" ]; then
    echo "Assembly failed or aborted (likely low reads). Creating empty bypass file."
    echo ">no_assembly_low_reads_myloasm" > "${SAMPLE_OUT}/assembly_primary.fa"
fi
echo "Finished: $(date)"
echo "Contigs: $(grep -c '>' ${SAMPLE_OUT}/assembly_primary.fa 2>/dev/null || echo 0)"
