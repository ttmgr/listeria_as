#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# -----------------------------------------------------------------------------
# Step 17 (v2 wrapper): run the v2 Python report generator.
# Input: compiled outputs under processing/
# Output: processing/report/pipeline_report.html
# Run: sbatch scripts/17_generate_report_v2.sh
# -----------------------------------------------------------------------------
WORK_DIR="/path/to/project"
echo "Generating pipeline report (v2)..."
echo "Start time: $(date)"
# Run the v2 script
python3 "${WORK_DIR}/scripts/17_generate_report_v2.py" "$WORK_DIR"
echo "Finished: $(date)"
echo "Report: ${WORK_DIR}/processing/report/pipeline_report.html"
