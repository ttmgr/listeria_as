#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=120G
#SBATCH -t 4:00:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=samtools_convert
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A_%a.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A_%a.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 1: BAM → FASTQ conversion using samtools
# Submit as: sbatch --array=1-66 01_samtools_bam2fastq.sh
# ============================================================

INPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/BAM"
OUTPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/samtools"
FILELIST="/lustre/groups/hpc/urban_lab/projects/tim/filelist.txt"

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
