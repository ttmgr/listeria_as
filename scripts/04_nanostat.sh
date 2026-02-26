#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=120G
#SBATCH -t 2:00:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=nanostat
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A_%a.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A_%a.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 4: Basic read statistics with NanoStat
# Submit as: sbatch --array=1-66 04_nanostat.sh
# ============================================================

INPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/nanofilt"
OUTPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/nanostat"
FILELIST="/lustre/groups/hpc/urban_lab/projects/tim/filelist.txt"

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
