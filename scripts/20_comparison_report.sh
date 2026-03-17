#!/bin/bash
# -----------------------------------------------------------------------------
# Step 20: Compare AS vs N results across the defined Black sample group.
# Input: summary tables from steps 15, 16, and 18
# Output: comparison tables and figures under processing/report
# Run: sbatch scripts/20_comparison_report.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

echo "Generating comparison report..."
echo "Start time: $(date)"
python3 "${SCRIPT_DIR}/20_comparison_report.py" "$WORK_DIR"
echo "Finished: $(date)"
