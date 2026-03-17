#!/bin/bash
# -----------------------------------------------------------------------------
# Step 7: Merge per-sample QC and Listeria summaries into shared tables.
# Input: processing/nanostat/*.txt and processing/listeria/listeria_summary.tsv
# Output: processing/stats/read_metrics_summary.csv and listeria_summary.csv
# Run: sbatch --dependency=afterok:<NANOSTAT_JOB_ID> scripts/07_compile_stats.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

NANOSTAT_DIR="${WORK_DIR}/processing/nanostat"
LISTERIA_DIR="${WORK_DIR}/processing/listeria"
OUTPUT_DIR="${WORK_DIR}/processing/stats"
mkdir -p "$OUTPUT_DIR"
echo "Compiling read statistics..."
echo "Start time: $(date)"
# ---- 1. Compile NanoStat results into CSV ----
OUTPUT_CSV="${OUTPUT_DIR}/read_metrics_summary.csv"
echo "sample,mean_read_length,mean_read_quality,median_read_length,median_read_quality,number_of_reads,read_length_N50,total_bases" > "$OUTPUT_CSV"
# Search for nanostat files with multiple naming patterns
shopt -s nullglob
STAT_FILES=("${NANOSTAT_DIR}"/nanostat_*.txt)
if [ ${#STAT_FILES[@]} -eq 0 ]; then
    echo "WARNING: No nanostat_*.txt files found. Trying *NanoStats.txt fallback..."
    STAT_FILES=("${NANOSTAT_DIR}"/*NanoStats.txt)
fi
if [ ${#STAT_FILES[@]} -eq 0 ]; then
    echo "WARNING: No nanostat_*.txt or *NanoStats.txt files. Trying bare files (NanoStat --tsv --name)..."
    for f in "${NANOSTAT_DIR}"/r*_barcode*; do
        [ -f "$f" ] && STAT_FILES+=("$f")
    done
fi
if [ ${#STAT_FILES[@]} -eq 0 ]; then
    echo "ERROR: No NanoStat output files found in ${NANOSTAT_DIR}"
    echo "  Tried: nanostat_*.txt, *NanoStats.txt, r*_barcode* (bare)"
    ls -la "${NANOSTAT_DIR}/" 2>/dev/null | head -20
fi
shopt -u nullglob

for stat_file in "${STAT_FILES[@]}"; do
    [ -f "$stat_file" ] || continue
    SAMPLE=$(basename "$stat_file")
    SAMPLE="${SAMPLE#nanostat_}"
    SAMPLE="${SAMPLE#filtered_}"
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
