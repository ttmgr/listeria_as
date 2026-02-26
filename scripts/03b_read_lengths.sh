#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=120G
#SBATCH -t 4:00:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=read_lengths
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A_%a.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A_%a.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 3b: Compute per-read lengths for raw and filtered FASTQ
# Produces aggregated TSVs consumed by addon step 26.
# Submit as: sbatch --array=1-66 --dependency=afterok:<NANOFILT_JOB> 03b_read_lengths.sh
# ============================================================

RAW_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/samtools"
FILT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/nanofilt"
WORK_DIR="/lustre/groups/hpc/urban_lab/projects/tim"
FILELIST="${WORK_DIR}/filelist.txt"

BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(basename "$BAM_FILE" .bam)

RAW_FASTQ="${RAW_DIR}/${BASENAME}.fastq"
FILT_FASTQ="${FILT_DIR}/filtered_${BASENAME}.fastq"

echo "Processing: ${BASENAME}"
echo "Start time: $(date)"

# Raw read lengths (aggregated: sample \t length \t state \t count)
if [ -f "$RAW_FASTQ" ]; then
    RAW_OUT="${WORK_DIR}/processing/read_lengths_raw_agg.tsv"
    TMP_RAW=$(mktemp)
    awk -v sample="$BASENAME" 'NR%4==2 {lens[length($0)]++} END {for (l in lens) print sample "\t" l "\t" "raw" "\t" lens[l]}' \
        "$RAW_FASTQ" > "$TMP_RAW"
    (
        flock -x 200
        cat "$TMP_RAW" >> "$RAW_OUT"
    ) 200>"${RAW_OUT}.lock"
    rm -f "$TMP_RAW"
    echo "Raw lengths written"
else
    echo "WARNING: Raw FASTQ not found: $RAW_FASTQ"
fi

# Filtered read lengths (aggregated)
if [ -f "$FILT_FASTQ" ]; then
    FILT_OUT="${WORK_DIR}/processing/read_lengths_filtered_agg.tsv"
    TMP_FILT=$(mktemp)
    awk -v sample="$BASENAME" 'NR%4==2 {lens[length($0)]++} END {for (l in lens) print sample "\t" l "\t" "filtered" "\t" lens[l]}' \
        "$FILT_FASTQ" > "$TMP_FILT"
    (
        flock -x 200
        cat "$TMP_FILT" >> "$FILT_OUT"
    ) 200>"${FILT_OUT}.lock"
    rm -f "$TMP_FILT"
    echo "Filtered lengths written"
else
    echo "WARNING: Filtered FASTQ not found: $FILT_FASTQ"
fi

echo "Finished: $(date)"
