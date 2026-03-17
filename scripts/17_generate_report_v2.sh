#!/bin/bash
# -----------------------------------------------------------------------------
# Step 17 (v2 wrapper): run the v2 Python report generator.
# Input: compiled outputs under processing/
# Output: processing/report/pipeline_report.html
# Run: sbatch scripts/17_generate_report_v2.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

echo "Generating pipeline report (v2)..."
echo "Start time: $(date)"
# Run the v2 script
python3 "${SCRIPT_DIR}/17_generate_report_v2.py" "$WORK_DIR"
echo "Finished: $(date)"
echo "Report: ${WORK_DIR}/processing/report/pipeline_report.html"
