#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=120G
#SBATCH -t 1:00:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=compile_stats
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 7: Compile all NanoStat results into a single CSV
# Submit as: sbatch --dependency=afterok:<NANOSTAT_JOB_ID> 07_compile_stats.sh
# ============================================================

NANOSTAT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/nanostat"
LISTERIA_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/listeria"
OUTPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/stats"

mkdir -p "$OUTPUT_DIR"

echo "Compiling read statistics..."
echo "Start time: $(date)"

# ---- 1. Compile NanoStat results into CSV ----
OUTPUT_CSV="${OUTPUT_DIR}/read_metrics_summary.csv"
echo "sample,mean_read_length,mean_read_quality,median_read_length,median_read_quality,number_of_reads,read_length_N50,total_bases" > "$OUTPUT_CSV"

for stat_file in "${NANOSTAT_DIR}"/barcode* "${NANOSTAT_DIR}"/nanostat_barcode*; do
    [ -f "$stat_file" ] || continue
    SAMPLE=$(basename "$stat_file")
    SAMPLE="${SAMPLE#nanostat_}"
    SAMPLE="${SAMPLE%.txt}"

    # Parse NanoStat TSV output (tab-separated: metric\tvalue)
    MEAN_LEN=$(grep "^mean_read_length" "$stat_file" | awk -F'\t' '{print $2}' || echo "NA")
    MEAN_QUAL=$(grep "^mean_qual" "$stat_file" | awk -F'\t' '{print $2}' || echo "NA")
    MED_LEN=$(grep "^median_read_length" "$stat_file" | awk -F'\t' '{print $2}' || echo "NA")
    MED_QUAL=$(grep "^median_qual" "$stat_file" | awk -F'\t' '{print $2}' || echo "NA")
    NUM_READS=$(grep "^number_of_reads" "$stat_file" | awk -F'\t' '{print $2}' || echo "NA")
    N50=$(grep "^n50" "$stat_file" | awk -F'\t' '{print $2}' || echo "NA")
    TOTAL_BASES=$(grep "^number_of_bases" "$stat_file" | awk -F'\t' '{print $2}' || echo "NA")

    echo "${SAMPLE},${MEAN_LEN},${MEAN_QUAL},${MED_LEN},${MED_QUAL},${NUM_READS},${N50},${TOTAL_BASES}" >> "$OUTPUT_CSV"
done

echo "Read metrics saved to: $OUTPUT_CSV"

# ---- 2. Sort and finalize Listeria summary ----
LISTERIA_CSV="${OUTPUT_DIR}/listeria_summary.csv"
echo "sample,listeria_reads,listeria_bases,listeria_mean_read_length,listeria_median_read_length" > "$LISTERIA_CSV"

if [ -f "${LISTERIA_DIR}/listeria_summary.tsv" ]; then
    # Keep the latest entry per sample in case summary.tsv was appended across reruns.
    awk -F'\t' '{last[$1]=$0} END {for (s in last) print last[s]}' "${LISTERIA_DIR}/listeria_summary.tsv" | \
        sort | tr '\t' ',' >> "$LISTERIA_CSV"
    echo "Listeria summary saved to: $LISTERIA_CSV"
else
    echo "WARNING: No listeria_summary.tsv found"
fi

# ---- 3. Print combined overview to stdout ----
echo ""
echo "========================================="
echo "         READ METRICS OVERVIEW           "
echo "========================================="
column -t -s',' "$OUTPUT_CSV"
echo ""
echo "========================================="
echo "         LISTERIA SUMMARY                "
echo "========================================="
if [ -f "$LISTERIA_CSV" ]; then
    column -t -s',' "$LISTERIA_CSV"
fi

echo ""
echo "Finished: $(date)"
