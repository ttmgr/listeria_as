#!/bin/bash
# -----------------------------------------------------------------------------
# Step 23: Generate cohort-level reports, split CSVs, and local plots.
# Runs all post-processing scripts that were previously local-only.
#
# Input: outputs from steps 15-22 (report CSVs, kraken2 CSVs, AMR overviews)
# Output:
#   - local_black_report_r{N}/ and local_blue_report_r{N}/ (cohort HTML reports)
#   - processing/kraken2_csv/{black,blue}/ (split Kraken CSVs)
#   - processing/amrfinder/overview/{black,blue}/ (split AMR CSVs)
#   - plots_black/ and plots_blue/ (local plots)
#
# Run: sbatch --dependency=afterok:<REPORT_JOB>:<KRAKEN_CSV_JOB> scripts/23_cohort_reports.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

echo "========================================="
echo "  Step 23: Cohort Reports & Post-Processing"
echo "========================================="
echo "Start time: $(date)"
echo ""

ERRORS=0

# ---- 1. Split Kraken2 CSVs by cohort ----
echo "--- Splitting Kraken2 CSVs by cohort ---"
if [ -f "${WORK_DIR}/processing/kraken2_csv/reads_classification.csv" ]; then
    python3 "${SCRIPT_DIR}/split_kraken_csv_by_cohort.py" "$WORK_DIR"
    [ $? -ne 0 ] && echo "WARNING: split_kraken_csv_by_cohort.py failed" && ((ERRORS++))
else
    echo "  Skipped: no reads_classification.csv yet"
fi
echo ""

# ---- 2. Split AMR overviews by cohort ----
echo "--- Splitting AMR overviews by cohort ---"
if [ -f "${WORK_DIR}/processing/amrfinder/overview/amr_reads_overview.csv" ]; then
    python3 "${SCRIPT_DIR}/split_amr_overview_by_cohort.py" "$WORK_DIR"
    [ $? -ne 0 ] && echo "WARNING: split_amr_overview_by_cohort.py failed" && ((ERRORS++))
else
    echo "  Skipped: no amr_reads_overview.csv yet"
fi
echo ""

# ---- 3. Cohort HTML reports (Black + Blue, per round) ----
echo "--- Generating cohort HTML reports ---"
python3 "${SCRIPT_DIR}/build_local_black_report.py" "$WORK_DIR"
[ $? -ne 0 ] && echo "WARNING: build_local_black_report.py failed" && ((ERRORS++))
echo ""

# ---- 4. Local plots (Black + Blue) ----
echo "--- Generating local plots ---"
python3 "${SCRIPT_DIR}/22_local_plots.py" "$WORK_DIR"
[ $? -ne 0 ] && echo "WARNING: 22_local_plots.py failed" && ((ERRORS++))
echo ""

# ---- Done ----
echo "========================================="
if [ "$ERRORS" -gt 0 ]; then
    echo "  Completed with $ERRORS warning(s)"
else
    echo "  All post-processing complete"
fi
echo "Finished: $(date)"
echo "========================================="

# List generated reports
echo ""
echo "Generated reports:"
ls -1 "${WORK_DIR}"/local_*_report*/report.html 2>/dev/null || echo "  (no cohort reports found)"
ls -1 "${WORK_DIR}"/plots_*/report.html 2>/dev/null || echo "  (no plot reports found)"
echo ""
echo "Pipeline reports:"
ls -1 "${WORK_DIR}/processing/report/"pipeline_report*.html 2>/dev/null
