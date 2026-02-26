#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=120G
#SBATCH -t 2:00:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=nanofilt
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A_%a.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A_%a.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 3: Length filtering with NanoFilt (remove reads < 100 bp)
# Submit as: sbatch --array=1-66 --dependency=afterok:<JOB2_ID> 03_nanofilt.sh
# ============================================================

INPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/porechop"
OUTPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/nanofilt"
FILELIST="/lustre/groups/hpc/urban_lab/projects/tim/filelist.txt"

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
