#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# ============================================================
# Step 8: metaMDBG assembly
# Submit as: sbatch --array=1-66 08_metamdbg.sh
# ============================================================
INPUT_DIR="/path/to/project/processing/nanofilt"
OUTPUT_DIR="/path/to/project/processing/mdbg"
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
metaMDBG asm \
    --out-dir "$SAMPLE_OUT" \
    --in-ont "$INPUT_FASTQ" \
    --threads 20
# Fallback: if metaMDBG aborted due to lack of reads, create a dummy file
if [ ! -f "${SAMPLE_OUT}/contigs.fasta.gz" ]; then
    echo "Assembly failed or aborted (likely low reads). Creating empty bypass file."
    echo ">no_assembly_low_reads_mdbg" > "${SAMPLE_OUT}/contigs.fasta"
    gzip "${SAMPLE_OUT}/contigs.fasta"
fi
echo "Finished: $(date)"
echo "Contigs: $(grep -c '>' ${SAMPLE_OUT}/contigs.fasta 2>/dev/null || echo 0)"
