#!/bin/bash
# -----------------------------------------------------------------------------
# Step 19: Build a fast reads-focused report before assembly-heavy steps finish.
# Input: read metrics, Listeria read summary, and AMR read summary
# Output: processing/report/reads_report.html
# Run: sbatch scripts/19_reads_report.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

echo "Generating reads-only report..."
echo "Start time: $(date)"
python3 "${SCRIPT_DIR}/19_reads_report.py" "$WORK_DIR"
echo "Finished: $(date)"
echo "Report: ${WORK_DIR}/processing/report/reads_report.html"
