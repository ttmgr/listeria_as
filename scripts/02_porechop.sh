#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=120G
#SBATCH -t 8:00:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=porechop
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A_%a.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A_%a.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 2: Adapter trimming with Porechop
# Submit as: sbatch --array=1-66 --dependency=afterok:<JOB1_ID> 02_porechop.sh
# ============================================================

INPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/samtools"
OUTPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/porechop"
FILELIST="/lustre/groups/hpc/urban_lab/projects/tim/filelist.txt"

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
