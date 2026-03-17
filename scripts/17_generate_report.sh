#!/bin/bash
# -----------------------------------------------------------------------------
# Step 17: Build the main HTML report for the full pipeline.
# Input: Listeria overview, AMR overview, and assembly stats outputs
# Output: processing/report/pipeline_report.html
# Run: sbatch --dependency=afterok:<OVERVIEW>:<AMR_OVERVIEW> scripts/17_generate_report.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

echo "Generating pipeline report..."
echo "Start time: $(date)"
python3 "${SCRIPT_DIR}/17_generate_report_v2.py" "$WORK_DIR"
echo "Finished: $(date)"
echo "Report: ${WORK_DIR}/processing/report/pipeline_report.html"
