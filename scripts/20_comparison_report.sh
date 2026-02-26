#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# -----------------------------------------------------------------------------
# Step 20: Compare AS vs N results across the defined Black sample group.
# Input: summary tables from steps 15, 16, and 18
# Output: comparison tables and figures under processing/report
# Run: sbatch scripts/20_comparison_report.sh
# -----------------------------------------------------------------------------
BASE_DIR="/path/to/project"
SCRIPT_DIR="$(dirname "$0")"
echo "Generating comparison report..."
echo "Start time: $(date)"
python3 "${SCRIPT_DIR}/20_comparison_report.py" "$BASE_DIR"
echo "Finished: $(date)"
