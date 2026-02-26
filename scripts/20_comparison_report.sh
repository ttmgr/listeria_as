#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# ============================================================
# Step 20: Black Sample Comparison Report (AS vs N)
# Depends on: Steps 15, 16, 18 (all summary data computed)
# ============================================================
BASE_DIR="/path/to/project"
SCRIPT_DIR="$(dirname "$0")"
echo "Generating comparison report..."
echo "Start time: $(date)"
python3 "${SCRIPT_DIR}/20_comparison_report.py" "$BASE_DIR"
echo "Finished: $(date)"
