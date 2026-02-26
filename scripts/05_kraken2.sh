#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=250G
#SBATCH -t 24:00:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=kraken2
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A_%a.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A_%a.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 5: Taxonomic classification with Kraken2
# Submit as: sbatch --array=1-66 05_kraken2.sh
# ============================================================


INPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/nanofilt"
OUTPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/kraken2"
KRAKEN2_DB="/lustre/groups/hpc/urban_lab/datasets/ncbi/kraken2_core"
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

KRAKEN2_REPORT="${OUTPUT_DIR}/report_${BASENAME}.txt"
KRAKEN2_OUTPUT="${OUTPUT_DIR}/classified_${BASENAME}.txt"

echo "Processing: ${BASENAME}"
echo "Input: ${INPUT_FASTQ}"
echo "Start time: $(date)"

kraken2 --db "$KRAKEN2_DB" \
        --use-names \
        --threads ${SLURM_CPUS_PER_TASK} \
        --report "$KRAKEN2_REPORT" \
        --output "$KRAKEN2_OUTPUT" \
        "$INPUT_FASTQ"

echo "Finished Kraken2: $(date)"
